---
title: "BigQuery Cloud Functions"
date: 2022-12-30T09:54:15-05:00
draft: false
---

One pattern that I learned about recently is that BigQuery allows for functions to be registered with itself that are backed by GCP Cloud Functions.  This is very powerful pattern that allows for much more complex functions than what you can do with the normal BigQuery User Defined Functions.  In this post I will walk through at a high level how to do with this Python but this can work with any language supported by Cloud Functions.


I found this example when attempting to call an rdkit function from BigQuery and discovered that someone had done this already in this [project](https://github.com/vjb-collab/cheminformatics-bq).  


The first step is to define the function.  For this example I will round trip a SMILES string from BigQuery to make sure it is a canonical SMILES using the rdkit library which also gives an example of third party dependency.


First lets write a function in a file called `function.py`.  The official (BigQuery documentation)[https://cloud.google.com/bigquery/docs/reference/standard-sql/remote-functions] has more details on this format but at a high level your python function receives a request object with a list of calls on it.  BigQuery will batch up rows together so you have to handle a list of calls.  Each call object has an array of values for each argument for a function.  Since this example only takes a single argumement we will just call the first one.  Finally you can call your function and return the values back as JSON in the format expected by BigQuery.


```python
import json
from rdkit import Chem


def rdkit_canonical_smiles(request):
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json["calls"]
        for call in calls:
            smiles = call[0]
            try:
                mol = Chem.MolFromSmiles(smiles)
                result = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                return_value.append(result)
            except Exception:
                return_value.append("")

        return_json = json.dumps({"replies": return_value}), 200
        return return_json
    except Exception:
        return json.dumps({"errorMessage": "something unexpected in input"}), 400
```

From there you need to create a requirements.txt file with the needed single dependency.

```text
rdkit==2022.3.5
```

Now you can register the fuction with cloud functions and BigQuery.  We will install it under the `cheminformatics` dataset that was created before.

```bash
SERVICE_ACCOUNT=$(bq show --location=US --format=prettyjson --connection "cheminformatics-connection" | jq -r '.cloudResource.serviceAccountId')

echo "Connection cheminformatics-connect service account: ${SERVICE_ACCOUNT}"

PROJ=$(gcloud config list --format 'value(core.project)')

PERM="roles/cloudfunctions.invoker"

TIMEOUT=600s
MEMORY=512MB
MAX_INSTANCES=100

gcloud beta functions deploy rdkit-canonical-smiles \
     --quiet --gen2 --region "us-east1" --entry-point rdkit_canonical_smiles --runtime python39 --trigger-http \
     --memory=$MEMORY --timeout=$TIMEOUT --max-instances=$MAX_INSTANCES  \
     --update-labels package=cheminformatics --update-labels function_type=remote_function --update-labels software_package=rdkit

CLOUD_TRIGGER_URL=$(gcloud beta functions describe rdkit-canonical-smiles --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

gcloud beta functions add-iam-policy-binding "rdkit-canonical-smiles" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} --gen2

gcloud run services add-iam-policy-binding "rdkit-canonical-smiles" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role="roles/run.invoker"

bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE or REPLACE FUNCTION cheminformatics.rdkit_canonical_smiles(smiles STRING) RETURNS STRING REMOTE WITH CONNECTION `us.cheminformatics-connection` OPTIONS (endpoint = @url, max_batching_rows = 2500)'
```




