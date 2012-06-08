#! /bin/bash
rm objs/libpbrt.a; make
pbrt ../../../hw3/hw3_bunnies.pbrt
exrtotiff hw3_bunnies.exr hw3_bunnies.tiff