#!/bin/sh

if [ ! -f libbito-core.so ] \
        || [ ! -f libbito.so ] \
        || [ ! -f doctest ] \
        || [ ! -f gp_doctest ]; then
    echo "Must run from cmake build directory after a successful build."
    exit 1
fi

ln -sf ../data .
ln -sf libbito.so bito.so
mkdir -p _ignore
export PYTHONPATH=.
pytest -s ../test/test_bito.py
./gp_doctest
./doctest
