#!/bin/bash
for i in {1..10..1}
do
   echo "processing:" $i
   ./build/run_DE_INT ../dataset/30/images_rectified $i output 5 1 12 12 0.05 0.95 0 8 0.40 0 1 0 2 0.25 2 1 0.05 0
done
