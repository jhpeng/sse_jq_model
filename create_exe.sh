#/bin/bash

mkdir include
mkdir lib

cd data_struct
make lib
cd ..
ar rcs lib/libdata.a data_struct/*.o
cp data_struct/*.h include

#rm -r include lib
