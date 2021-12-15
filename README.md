# BDIS: Bayesian Dense Inverse Searching Methodfor Real-Time Stereo Surgical Image Matching #

This is the CPU level real-time stereo software for minimally invasive surgery.

 
  
## Compiling ##

The program was only tested under a 64-bit Linux distribution.
SSE instructions from built-in X86 functions for GNU GCC were used.


```
cd OF_dis_bayesian
mkdir build
cd build
cmake ../
make -j
```

The code depends on Eigen3 and OpenCV.
      

## Usage ##
1. ./run_bash_*.sh generates the disparity.   (Use 'chmod +x run_bash_*.sh' to assign permission )   
2. Run Main_batch.m (data set 6, 20 and 21) or Main_synthetic.m (data set 30, 31, 32 and 33) in the matlab folder for visualization.      
3. [Optional]. To test the coverage rate of the variance, please modify the last parameter in bash_*.sh from 0 to 1.       
      



## LICENCE CONDITIONS ##

This work is released under GPLv3 license. For commercial purposes, please contact the authors: jingweisong@yahoo.com











