# Bayesian dense inverse searching algorithm for real-time stereo matching in minimally invasive surgery #

This is the CPU level real-time stereo software for minimally invasive surgery.

 
  
## Compiling ##

The program was only tested under a 64-bit Linux distribution.
SSE instructions from built-in X86 functions for GNU GCC were used.


```
Compile:
chmod +x build.sh
./build.sh
```

The code depends on Eigen3 and OpenCV.
      

## Usage ##
For testing:      
1. ./run_bash_*.sh generates the disparity.      
2. Run Main.m in the matlab folder for visualization.    

To integrate into your project, just complile include 'getDisp.h' and 'libFastPatchOF.so'  in your project.
      

Citation:      
@article{song2022bdis,      
  title={{BDIS}: Bayesian Dense Inverse Searching Method for Real-Time Stereo Surgical Image Matching},      
  author={Song, Jingwei and Zhu, Qiuchen and Lin, Jianyu and Ghaffari, Maani},      
  journal={IEEE Transactions on Robotics},      
  volume={39},      
  number={2},      
  pages={1388--1406},      
  year={2022},      
  publisher={IEEE}      
}      
@inproceedings{song2022bayesian,      
  title={Bayesian dense inverse searching algorithm for real-time stereo matching in minimally invasive surgery},      
  author={Song, Jingwei and Zhu, Qiuchen and Lin, Jianyu and Ghaffari, Maani},      
  booktitle={International Conference on Medical Image Computing and Computer-Assisted Intervention},      
  pages={333--344},      
  year={2022},      
  organization={Springer}      
}      


## LICENCE CONDITIONS ##

This work is released under GPLv3 license. For commercial purposes, please contact the authors: jingweisong@yahoo.com, jingweisong.eng@outlook.com










