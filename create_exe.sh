#/bin/bash

mkdir include
mkdir lib

cd data_struct
rm *.o
make lib
cd ../
ar rcs lib/libds.a data_struct/*.o
cp data_struct/*.h include

cd algorithm
rm *.o
make lib
cd ../
ar rcs lib/libal.a algorithm/*.o
cp algorithm/*.h include


#rm -r include lib
