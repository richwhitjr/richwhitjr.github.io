---
title: "Bazel + Python + CUDA"
date: 2022-12-27T09:13:30-05:00
draft: false
---

One of the gaps with the current [Bazel Python Rules](https://github.com/bazelbuild/rules_python) is that it does
not handle pip dependencies well when you need to install a different library if you have cuda available or not. 
This is a problem with torch where a different install is needed depending on cuda.  This post is a fairly minimal
workaround for the problem where a new repository rule is created that uses different requirements.txt files if CUDA
is found or not.  It uses the new [Cuda Rules](https://github.com/bazel-contrib/rules_cuda) to do this.


First in your `WORKSPACE` make sure that python is installed and you pick the correct toolchain for your version.  Here I
am using python 3.10 and the version of the toolchain that was current for when this was posted.

```python
rules_python_version = "834149dfcd9e0dcb9d713caeb6bf5b0584601392"

http_archive(
    name = "rules_python",
    sha256 = "0df35e6e46ea7abe0ddacc32a1d2f8884fc13ef00654411ca97e9f2c9b0589dc",
    strip_prefix = "rules_python-{}".format(rules_python_version),
    url = "https://github.com/bazelbuild/rules_python/archive/{}.zip".format(rules_python_version),
)

load("@rules_python//python:repositories.bzl", "python_register_toolchains")

python_register_toolchains(
    name = "python3_10",
    python_version = "3.10",
)

load("@python3_10//:defs.bzl", "interpreter")
```

Next in the same `WORKSPACE` file also install the rules cuda.  These rules will be used to check if cuda is installed on the machine
running the workspace.

```python
http_archive(
    name = "rules_cuda",
    sha256 = "d1ec8d5cbec3514c6ed00eabf168a4b4e0b9dc07797aea98671722bc257dcc23",
    strip_prefix = "rules_cuda-0c32b04d8af865182a3b6367c1bdd051a699c7de",
    urls = ["https://github.com/bazel-contrib/rules_cuda/archive/0c32b04d8af865182a3b6367c1bdd051a699c7de.tar.gz"],
)

load("@rules_cuda//cuda:repositories.bzl", "register_detected_cuda_toolchains", "rules_cuda_dependencies")

rules_cuda_dependencies()

register_detected_cuda_toolchains()
```

After this we need to create a new rule that allows for the swapping different requirements file based on CUDA availablility.
In this example the rule will be created in at `tools/python/defs.bzl`.


At a high level the rule works like this.  Take two different requirements files and call `pip_install` on them.  This will install
two different requirements like normal using the python rules.  To do the switch however the new rule generates a new pip rules that refers
to the correct requirements file based on cuda being available or not.  The complete code is seen below:

```python
load("@rules_python//python:pip.bzl", "pip_install")
load("@rules_cuda//cuda/private:repositories.bzl", "detect_cuda_toolkit")

def _multi_pip_parse_impl(rctx):
    """Rule that finds the correct requirements file if cuda is available"""
    cuda = detect_cuda_toolkit(rctx)
    deps = "cpu_deps"
    if cuda.path:
        deps = "cuda_deps"
    requirements_bzl = """
def _clean_name(name):
    return name.replace("-", "_").replace(".", "_").lower()
def requirement(name):
    return "@{name}//pypi__" + _clean_name(name)
""".format(
        name = deps,
    )
    rctx.file("requirements.bzl", requirements_bzl)
    rctx.file("BUILD.bazel", "exports_files(['requirements.bzl'])")

_multi_pip_parse = repository_rule(_multi_pip_parse_impl)

def multi_pip(name, python_interpreter_target, cuda, cpu):
    """Install multiple different requirements files based on system properties"""

    # Install CPU Deps
    pip_install(
        name = "cpu_deps",
        timeout = 60 * 30,  # 30 minutes
        python_interpreter_target = python_interpreter_target,
        requirements = cpu,
    )

    # Install CUDA Deps
    pip_install(
        name = "cuda_deps",
        timeout = 60 * 30,  # 30 minutes
        python_interpreter_target = python_interpreter_target,
        requirements = cuda,
    )

    # Pick the correct dependency based on finding CUDA
    return _multi_pip_parse(name = name)
```

Finally we need to make use of new rule back in the `WORKSPACE` file.  This will load 
the new rule and give it two different requirements files wit the interpreter defined above.

```python
load("//tools/python:defs.bzl", "multi_pip")

multi_pip(
    name = "deps",
    cpu = "//build:cpu_requirements.txt",
    cuda = "//build:gpu_requirements.txt",
    python_interpreter_target = interpreter,
)
```

Now we can use this like normal python rules to pull in pip dependencies using the name of given
to the `multi-pip` library:

```python
load("@rules_python//python:defs.bzl", "py_binary")
load("@deps//:requirements.bzl", "requirement")

py_binary(
    name = "main",
    srcs = ["main.py"],
    deps = [requirement("torch")],
)
```