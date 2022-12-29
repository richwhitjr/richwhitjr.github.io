---
title: "Bazel + Beam/Dataflow"
date: 2022-12-27T20:01:55-05:00
draft: false
---

Python Beam on DataFlow allows for a [few options](https://beam.apache.org/documentation/sdks/python-pipeline-dependencies/) for including internal and external dependencies in your jobs.  None of these works very well with bazel however because of the way third party dependencies are packaged.


The main issue is that Bazel Python does not install pip libraries on the system or in the container if you choose that route.  That means there is not a great
way to register external dependencies with the job like mentioned in the above document.  The one route that I have found to work well and that is described below is to instead use a custom image through Bazel but with a custom launch script that sets up the pip dependencies correctly on the `PYTHONPATH`.


The first step is to create the start script called by the Docker image.  We will create this in `start.sh`.  This script is doing the work on walking through the pip libraries as installed by [Bazel Docker](https://github.com/bazelbuild/rules_docker#py3_image).

```bash
SCRIPT_PATH=$PWD

for f in `ls -d external/*deps/pypi__*/site-packages`; do
    if [ -z "$PYTHONPATH" ]
    then
        export PYTHONPATH=$SCRIPT_PATH/$f
    else
        export PYTHONPATH=$PYTHONPATH:$SCRIPT_PATH/$f
    fi    
done;

PYTHONPATH=$PYTHONPATH:$SCRIPT_PATH /opt/apache/beam/boot "$@"
```

The next steps are following the Beam documentation where a Dockerfile is created with a custom entry script.  We copy the script we created above:

```text
FROM apache/beam_python3.10_sdk:2.43.0

COPY start.sh /opt/start.sh

RUN chmod 777 /opt/start.sh

ENTRYPOINT ["/opt/start.sh"]
```

After the image is pushed to somewhere like dockerhub we need register it in the Bazel `WORKSPACE` file like any other docker image pull.  See the [Github page](https://github.com/bazelbuild/rules_docker/blob/master/docs/container.md#container_pull) for more details.

```python
container_pull(
    name = "dataflow_container",
    registry = "....",
    repository = "....",
    tag = "latest",
)
```

Finally we can create push the create image to a Docker repository using the Bazel rules `bazel run //dataflow:push` if the follow build is located at `dataflow/BUILD.bazel`.

```python
load("@rules_python//python:defs.bzl", "py_binary")
load("@io_bazel_rules_docker//python3:image.bzl", "py3_image")
load("@io_bazel_rules_docker//container:container.bzl", "container_image")
load("@deps//:requirements.bzl", "requirement")

py_binary(
    name = "main",
    srcs = ["main.py"],
    deps = [requirement("torch")]
)

py3_image(
    name = "py_image",
    srcs = ["main.py"],
    main = "main,
    base = "@dataflow_container//image",
    deps = [requirement("torch")]
)

container_image(
    name = "image",
    base = "py_image",
    entrypoint = ["/opt/start.sh"],
)

container_push(
    name = "push,
    format = "Docker",
    image = ":image",
    registry = "...",
    repository = "...,
    tag = "latest",
)
```

Now to run the job we need execute the main bazel target but specifying the docker image we just pushed pointing the job at the Docker repository.

```bash
bazel run //dataflow:push
bazel run //dataflow:main -- \
  --runner DataflowRunner \
  --sdk_container_image=$REPO/$IMAGE:latest
```