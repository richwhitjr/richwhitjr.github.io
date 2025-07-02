---
title: "RDkit + GCP Sharding"
date: 2025-07-02T11:18:41-04:00
draft: false
---

## Sharding Your RDKit-Backed Chemical Database on GCP

*Distributed Chemical Informatics | RDKit & Google Cloud Run*

---

### Introduction

When you’re working with large libraries of molecules (e.g., millions of compounds encoded as SMILES or SDF), in-memory lookups can become a bottleneck. RDKit makes substructure, fingerprint, and similarity searches effortless—but how do you scale that across many machines without reinventing the wheel?

In this post, we’ll:

1. **Pre-shard** an RDKit-compatible chemical dataset
2. Build a **per-shard microservice** in Cloud Run
3. Deploy **N independent shard services**
4. Implement a **fan-out aggregator**
5. (Optionally) front it all with **API Gateway**

---

## 1. Pre-shard Your Dataset

Splitting your SDF/SMILES file into *N* self-contained shards lets each service load only its slice:

```bash
# Example: split by line count (for SMILES+ID file)
split -l $((TOTAL_LINES/N)) all_molecules.smi shard-
# Rename:
for i in shard-*; do mv "$i" "shard-${i#shard-}.smi"; done
# Upload:
gsutil cp shard-*.smi gs://my-chem-db/
```

> **Pro tip:** If you have SDFs, you can use RDKit’s Python API (`Chem.SDMolSupplier`) to write out evenly sized chunks.

---

## 2. Build the Shard Microservice

Each shard service loads **exactly one** `.smi` or `.sdf` at startup and uses RDKit for queries:

```python
# app.py
import os, json, flask
from rdkit import Chem
from rdkit.Chem import AllChem

app = flask.Flask(__name__)
SHARD_ID    = int(os.getenv("SHARD_ID"))
SHARD_COUNT = int(os.getenv("SHARD_COUNT"))
BUCKET      = os.getenv("SHARD_BUCKET")

# 1. Download your shard file on cold start
from google.cloud import storage
client = storage.Client()
blob = client.bucket(BUCKET).blob(f"shard-{SHARD_ID}.smi")
blob.download_to_filename("/tmp/shard.smi")

# 2. Load molecules into memory
mols = []
with open("/tmp/shard.smi") as f:
    for line in f:
        smi, cid = line.strip().split()
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol.SetProp("_Name", cid)
            AllChem.Compute2DCoords(mol)
            mols.append(mol)

@app.route("/lookup")
def lookup():
    query_smiles = flask.request.args["q"]
    qmol = Chem.MolFromSmiles(query_smiles)
    fps = AllChem.GetMorganFingerprintAsBitVect(qmol, 2)
    results = []
    for m in mols:
        score = DataStructs.TanimotoSimilarity(
            fps,
            AllChem.GetMorganFingerprintAsBitVect(m, 2),
        )
        if score > 0.7:
            results.append({"id": m.GetProp("_Name"), "score": score})
    return flask.jsonify(results)

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080)
```

* **Env vars**: `SHARD_ID`, `SHARD_COUNT`, `SHARD_BUCKET`
* **RDKit**: substructure or fingerprint search out of the box

---

## 3. Deploy N Shard Services to Cloud Run

Use the same container image, but parameterize each deployment:

```bash
for i in $(seq 0 $((N-1))); do
  gcloud run deploy chem-shard-$i \
    --image gcr.io/PROJECT/chem-shard:latest \
    --region us-central1 \
    --set-env-vars SHARD_ID=$i,SHARD_COUNT=$N,SHARD_BUCKET=my-chem-db \
    --allow-unauthenticated
done
```

* **Auto-scaling**: each shard scales to zero when idle
* **Caching**: the in-memory list `mols` stays hot

> *Tip:* if you want a minimum of M total instances across all shards, set `--min-instances=$((M/N))`.

---

## 4. Write the Aggregator Service

This Cloud Run service fans out incoming queries to all N shards **in parallel**, merges the results, and returns them:

```python
# aggregator.py
import os, asyncio, aiohttp, flask

app = flask.Flask(__name__)
SHARD_COUNT = int(os.getenv("SHARD_COUNT"))
SHARD_URLS  = [
    os.getenv(f"SHARD_URL_{i}") for i in range(SHARD_COUNT)
]

async def query_shard(session, url, q):
    async with session.get(url, params={"q": q}) as resp:
        return await resp.json()

@app.route("/chem")
def chem_lookup():
    q = flask.request.args["q"]
    results = asyncio.run(_fan_out(q))
    # Flatten, sort by score descending, dedupe if needed
    merged = sorted(
        [item for sub in results for item in sub],
        key=lambda x: x["score"], reverse=True
    )
    return flask.jsonify(merged)

async def _fan_out(q):
    async with aiohttp.ClientSession() as session:
        tasks = [
            query_shard(session, url + "/lookup", q)
            for url in SHARD_URLS
        ]
        return await asyncio.gather(*tasks)

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080)
```

Deploy it similarly, passing each shard’s URL:

```bash
gcloud run deploy chem-aggregator \
  --image gcr.io/PROJECT/chem-aggregator:latest \
  --set-env-vars SHARD_COUNT=$N,\
SHARD_URL_0=https://chem-shard-0-xxxx.a.run.app,\
SHARD_URL_1=https://chem-shard-1-xxxx.a.run.app,...
```

---

## 5. (Optional) Single Endpoint via API Gateway

For a unified hostname, auth, and rate limiting:

1. **Define** an OpenAPI spec that routes `/chem` to your aggregator.
2. **Create** an API config and gateway:

   ```bash
   gcloud api-gateway api-configs create chem-config \
     --api=chem-api --openapi-spec=openapi.yaml \
     --project=PROJECT

   gcloud api-gateway gateways create chem-gateway \
     --api=chem-api --api-config=chem-config \
     --location=us-central1
   ```

Clients can now hit `https://chem-gateway-…/chem?q=<SMILES>`.

---

## Conclusion & Next Steps

By combining RDKit’s power with GCP’s serverless footprint, you get:

* **High throughput**: in-memory shards and concurrent fan-out
* **Cost-efficiency**: pay-per-use, scale-to-zero
* **Simplicity**: a handful of services, clear separation of concerns

**Future enhancements** might include:

* **Cloud Workflows** for orchestration
* **Pub/Sub + Cloud Functions** for async, event-driven queries
* A migration onto **Bigtable** or **Spanner** for fully managed sharding
