Installation Guide
==================

## C++ Interface

### Compile C++ sources and run tests
---

To compile the tests of the C++ interface simply

```bash
mkdir -p test/build && cd test/build
cmake ..
make
```

In WSL (Windows Subsystem Linux), you can run the following command to install libc6-dev-i386. This will be required for `ieeefp.h` which is used by `qd` library,

    sudo apt-get install libc6-dev-i386

Also to install `mkl` related dependencies, run the following,

    wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
    sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
    sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
    sudo apt-get update
    sudo apt-get install intel-mkl-2020.4-912
    sudo sh -c "echo /opt/intel/mkl/lib/intel64 > /etc/ld.so.conf.d/intel-mkl.conf"
    sudo ldconfig
    export CPLUS_INCLUDE_PATH="/opt/intel/mkl/include:$CPLUS_INCLUDE_PATH"

You can run the tests by `cmake test` or `ctest -jK` where `K` the number of `CPU` threads. By adding the option `--verbose` to `ctest` you get more information about the tests, *e.g.* time per test, volume computed and the name of the polytope or convex body.

### Development environment from Docker container
---
Optionally, it is possible to setup a docker contaner with development environment. To get started with docker, see
[here](https://docs.docker.com/get-started/). Below, here is how a Dockerfile can be written
```
FROM ubuntu:18.04

RUN apt-get update && apt-get install -y g++ cmake lp-solve && \
    rm -rf /var/lib/apt/lists/*
```

The user should create a file `Dockerfile.dev` with above content inside a temporary folder (e.g. `docker`) then create an image and run a container:

```bash
# build docker image
cd docker
docker build -t volesti:dev -f Dockerfile.dev .
# check built image
docker images | grep volesti
# run a container in an interactive mode from volesti source folder
docker run -it -v $PWD:/volesti -w /volesti --name=volesti-dev volesti:dev /bin/bash
```

## R Interface

An ``R`` interface is available from the package [Rvolesti](https://github.com/GeomScale/Rvolesti).

## Python Interface

A ``python`` interface is available from the package [dingo](https://github.com/GeomScale/dingo).


