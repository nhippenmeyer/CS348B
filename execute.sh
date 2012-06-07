#! /bin/bash

rm lightfield.txt
rm lightfield.bin
rm objs/libpbrt.a; make
pbrt ../../../hw3/hw3_bunnies.pbrt --ncores 1
exrtotiff hw3_bunnies.exr hw3_bunnies.tiff