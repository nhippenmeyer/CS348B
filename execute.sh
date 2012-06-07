#! /bin/bash

rm lightfield.txt
rm lightfield.bin
rm objs/libpbrt.a; make
time pbrt ../../../hw3/hw3_bunnies.pbrt
exrtotiff hw3_bunnies.exr hw3_bunnies.tiff