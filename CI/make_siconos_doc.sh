#!bin/bash
apt-get update -qq && apt-get install -y -qq cmake git-core wget make \
			      cmake libboost-dev libgmp-dev swig gcc gfortran g++ liblapack-dev libatlas-base-dev \
			      lp-solve liblpsolve55-dev python3-dev libpython3-dev bash swig doxygen python3-dev python3-pip

pip3 install -r ./docs/requirements.txt
mkdir build
cd build
cmake ../ -DUSER_OPTIONS_FILE=$HOME/siconos-mirror/CI/siconos_docs.cmake
make doxygen
make doc
