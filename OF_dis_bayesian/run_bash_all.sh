# dataset 6
#!/bin/bash
for i in {1..1201..50}
do
   echo "processing:" $i
   ./build/run_DE_INT ../dataset/6/images_rectified $i output 5 1 12 12 0.05 0.95 0 10 0.55 0 1 0 2 0.75 2 0.5 0.05 0

done
# dataset 20
#!/bin/bash
for i in {801..2401..100}
do
   echo "processing:" $i
   ./build/run_DE_INT ../dataset/20/images_rectified $i output 5 1 12 12 0.05 0.95 0 10 0.60 0 1 0 2 0.25 2 0.5 0.05 0
done
# dataset 21
#!/bin/bash
for i in {1001..4201..100}
do
   echo "processing:" $i
   ./build/run_DE_INT ../dataset/21/images_rectified $i output 5 1 12 12 0.05 0.95 0 10 0.60 0 1 0 2 0.25 2 0.5 0.05 0
done
# dataset 30
#!/bin/bash
for i in {1..10..1}
do
   echo "processing:" $i
   ./build/run_DE_INT ../dataset/30/images_rectified $i output 5 1 12 12 0.05 0.95 0 8 0.40 0 1 0 2 0.25 2 1 0.05 0
done
# dataset 31
#!/bin/bash
for i in {1..16..1}
do
   echo "processing:" $i
   ./build/run_DE_INT ../dataset/31/images_rectified $i output 5 1 12 12 0.05 0.95 0 8 0.40 0 1 0 2 0.25 2 1 0.05 0
done
# dataset 32

for i in {1..26..1}
do
   echo "processing:" $i
   ./build/run_DE_INT ../dataset/32/images_rectified $i output 5 1 12 12 0.05 0.95 0 8 0.40 0 1 0 2 0.25 2 1 0.05 0
done
# dataset 33
#!/bin/bash
for i in {1..26..1}
do
   echo "processing:" $i
   ./build/run_DE_INT ../dataset/33/images_rectified $i output 5 1 12 12 0.05 0.95 0 8 0.40 0 1 0 2 0.25 2 1 0.05 0
done
