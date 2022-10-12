# BDIS: Bayesian Dense Inverse Searching Methodfor Real-Time Stereo Surgical Image Matching #

This is the CPU level real-time stereo software for minimally invasive surgery. The paper has been accepted by IEEE Transaction on Robotics (https://arxiv.org/abs/2205.03133).  

 
  
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

Parameter setting:    

1.Coarsest scale                               (here: 5)    
2.Finest scale                                 (here: 3)    
3/4.Min./Max. iterations                       (here: 12)    
5./6./7. Early stopping parameters    
8.Patch size                                   (here: 8)    
9.Patch overlap                                (here: 0.4)    
10.Use forward-backward consistency             (here: 0/no)    
11.Mean-normalize patches                       (here: 1/yes)    
12.Cost function                                (here: 0/L2)  Alternatives: 1/L1, 2/Huber, 10/NCC    
13.Verbosity                                   (here: 2) Alternatives: 0/no output, 1/only flow runtime, 2/total runtime    
14.ratio_patch_valid:                       (here: 0.75) Minimal ratio of valid patch. (ratio of the valid points in the patch)     
15.num_window:                               (here: 2) number of Bayesian window (left)    
16.unit_disturb:                               (here: 0.5) Gap between consectutivePSe window    
17.min_prob:                                  (here: 0.05) Minimum probability for the output    
18.bool_var:                                    (here: False) Export variance    
      



## LICENCE CONDITIONS ##

This work is released under GPLv3 license. For commercial purposes, please contact the authors: jingweisong@yahoo.com











