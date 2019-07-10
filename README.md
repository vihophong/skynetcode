# skynetcode
An user-defined code written in C++ based on SkyNet framework

SkyNet is available at: https://bitbucket.org/jlippuner/skynet/

Lippuner, Jonas and Roberts, Luke F. (2017) SkyNet: A Modular Nuclear Reaction Network Library. Astrophysical Journal Supplement Series, 233 (2)

## Installation:

Requirement: cmake3 and skynet code installed 

Define enviroment variables in ~/.bashrc

export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64:$LD_LIBRARY_PATH

export PYTHONPATH=/home/phong/projects/SkyNet/install/skynet-install/lib:$PYTHONPATH

export SKYNETSYS=/home/phong/projects/SkyNet/install/skynet-install/lib/cmake/SkyNet

export MKLROOT=/opt/intel/compilers_and_libraries/linux/mkl


Go to parent directory


mkdir build-skynetcode

cd build-skynetcode

cmake ../skynetcode

make
