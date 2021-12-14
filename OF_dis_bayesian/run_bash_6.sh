#!/bin/bash
for i in {1..1201..50}
do
   echo "processing:" $i
   ./build/run_DE_INT ../dataset/6/images_rectified $i output 5 1 12 12 0.05 0.95 0 10 0.55 0 1 0 2 0.75 2 0.5 0.05 0

done
