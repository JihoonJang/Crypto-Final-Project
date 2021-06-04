# Cryptography

```bash
git submodule update --init --recursive

# Install HElib
cd HElib
mkdir build
cd build
cmake -DPACKAGE_BUILD=ON -DCMAKE_INSTALL_PREFIX=../../ ..
make -j16
make install

# Run program
cd ../../matrix_multiplication/
mkdir build
cd build
cmake -Dhelib_DIR=../../helib_pack/share/cmake/helib ../ ..
make
./bin/matmul
```