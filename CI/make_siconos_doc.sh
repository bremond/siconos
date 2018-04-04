#!bin/bash
pip install -r ./docs/requirements.txt
mkdir build
cd build
cmake ../ -DUSER_OPTIONS_FILE=../CI/siconos_docs.cmake
make doxygen
make doc
