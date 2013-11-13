#! /bin/bash
#Compiles DFM's Python implementation of the BLS library.
#Uses the inplace operation to circumvent the need for root.

echo -e 'Removing the previous builds...\n\n\n'
rm -rf ./build
rm -rf bls/_*.so

echo -e 'Setting up new build.\n\n\n'
python setup.py build_ext --inplace

echo -e 'Done.'
