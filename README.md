NOTE: The program is still under development, and some bugs still exist in PPPLib. Comments and suggestions from users 
      are welcome, please send them to the author.

**INTRODUCTION**

Precise Point Positioning Library (PPPLib) is a multi-GNSS data processing software designed to
process multi-frequency data from GPS, BeiDou, Galileo, GLONASS and QZSS. PPPLib is written
in the C/C++ programming language. It can compile, run on both Linux and Windows operating
system. PPPLib mainly performs precise point positioning from single- to triple-frequency based
on ionosphere-free or uncombined observations. Moreover, it provides abundant solutions including
positioning, tropospheric delay, ionospheric delay as well as ambiguity information. Useful scripts
and visualization tool are also provided for data download, batch processing or solution presentation.
We give a preliminary review, including positioning accuracy and convergence time of PPP using
dual-frequency, ionosphere-free from single-system to multi-GNSS, to show the working status of
current version of the software. In addition, the software also supports post-processing kinematic
mode and INS/GNSS loosely coupled mode for real kinematic scene. Detailed information can be found in User Manual.

**FEATURES**

* Supports multi-GNSS or their combinations with each other
* Standard single point positioning
* Single-frequency ionosphere-free PPP
* Dual-frequency ionosphere-free/uncombined PPP
* Triple-frequency ionosphere-free with one combination or two combination/uncombined PPP
* Multi-GNSS post-kinematic processing
* INS/GNSS loosely coupled
* Support BDS-3 new satellites and signals
* Convenient visualization

**WILL BE**　

1. GNSS/INS tightly coupled
2. PPP ambiguity resolution
3. Robust processing strategy for real-kinematic scene
4. Enhance algorithms for GNSS/INS integration

**HOW TO USE**
  
  PPPLib uses CMAKE for project management. Currently the software supports compile and run in
  Linux and Windows. We have tested the code under the virtual machine Ubantu 16.0, Ubantu 16.0
  and Win10. PPPLib uses the [easyloggingpp](https://github.com/amrayn/easyloggingpp) library to log information and [Eigen](https://github.com/artsy/eigen) is used to perform
  matrix operations. Both of them are single header-only library for C++ applications. If there are any
  problems about easyloggingpp and Eigen library, please refer to the website
  https://github.com/amrayn/easyloggingpp and https://github.com/artsy/eigen for help. Note that
  the author uses JetBrains CLion2019.3 software platform for developing and debug under
  Ubantu16.0. The CMAKE version is 3.15.3 and the C++ complier is gcc 5.4.0. 
   
* Linux user
    
  You can build and use the project follow steps 0-6:
  
  1. open Terminal, and cd [your path]
  2. git clone https://github.com/heiwa0519/PPPLib.git
  3. cd PPPLib
  4. mkdir build and cd [your path]/PPPLib/build
  5. cmake ..
  6. make
 
  After four steps, the executable files are generated in the [your path]/PPPLib/bin folder. You can use
  PPPLib in Terminal like:
  
  1. cd ../bin
  2. ./PPPLibMain -C ../conf/PPPLib.ini -M PPP-KINE -S G -L 128 -T 2019/12/01
  
  Note that the path in configuration file should set to your local dir.  

* Windows user

  You should install the MinGW on your Windows computer for gcc compiler. It is recommended to
  use CLion IDE for run and debug. Load the PPPLib project, after cmake completed, you can add program arguments “-C ../conf/PPPLib.ini
  -M PPP-KINE -S G -L 128 -T 2019/12/01” in Run/Debug Configurations and then the project is ready for running.
  
**DATASET**
  
  PPPLib provides different dataset to evaluate its performance. This dataset contains three kinds of
  data, namely static GNSS data, suburban vehicle data and urban data, which can be used to evaluate
  the performance of PPP, PPK, GNSS/INS integration mode of the software. Note some data are
  collected from the Internet, if you use the dataset, please cite related paper or indicate the source.
  For mainland China users, please download the dataset using the [Baidu Cloud link(8888)](https://pan.baidu.com/s/1Yr5g9O_U52Wp7T_K2aPoqg#list/path=%2F).
  Also can be found in repository [PPPLib-Dataset](https://github.com/heiwa0519/PPPLib-Dataset).
  
**ACKNOWLEDGEMENT**

 First of all, I pay tribute to Mr. Tomoji Tkasu, the author of [RTKLIB](https://github.com/tomojitakasu/RTKLIB) software. I admire him for his 
 selfless open source spirit and elegant programming. The most functions of PPPLib are refer to
 RTKLIB. The software is also refers to [rtkexplorer](http://rtkexplorer.com/), HPRTK of Mowen Li from Shandong university,
 the carvig software of Jinlan Su from Wuhan University, the GAMP software of Feng Zhou from
 Shangdong University of Science and Technology, the Multi-Sensor Fusion software of Wentao
 Fang from Wuhan University, the MG-APP software of Gongwei Xiao from Institute of Geodesy
 and Geophysics, Chinese Academy of Sciences and the PINS software of Gongmin Yan from
 Northwestern Polytechnic University. The easyloggingpp is used to log information and Eigen is
 used to perform matrix operations. Thanks to the authors of above software. Many thanks are due
 to [Steve Hillia](https://www.researchgate.net/profile/Steve_Hilla) from Notional Oceanic and Atmospheric Administration for his detailed suggestions.
  
**CONTACT AUTHOR**
  
  Any suggestions, corrections, and comments about PPPLib are sincerely welcomed and could be
  sent to:
  * Author: Chao Chen
  * QQ: 565681993
  * E-mail: cchen@cumt.edu.cn
  * Addres: School of Environment and Geo-informatics, China University of Mining and Technology



**SUPPORTED BY [Guobin Chang's Lab](https://www.researchgate.net/lab/Guobin-Chang-Lab)**