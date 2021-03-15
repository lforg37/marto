# Marto: Modern ARithmetic Tools for high-level synthesis

Marto provides C++ headers to implement custom sized arithmetic objects such as:

+ Custom sized posits and their environment (including the quire)
+ Custom sized IEEE-754 numbers 
+ Custom sized Kulisch accumulators (and sums of products)


# Marto usage

Marto is a header only library.
It is based on [Hint](https://github.com/yuguen/hint) in order to provide support for many HLS backends.

# Building unit tests 

+ Make sure to have the [Hint](https://github.com/yuguen/hint) library installed 
+ Make sure to have the [Softposit](https://gitlab.com/cerlane/SoftPosit) library installed
+ Make sure to have the ```ap_int``` library installed
+ Clone the current repository 
+ Create a build directory 
+ Create a cmake build 
+ Make the tests
+ Launch the tests

```Shell
git clone https://gitlab.inria.fr/lforget/marto 
cd marto
mkdir build
cd build 
cmake .. -DSOFTPOSIT_H=<PATH_TO_SOFTPOSIT_INSTALL_DIR>/source/include
         -DSOFTPOSIT_LIB=<PATH_TO_SOFTPOSIT_INSTALL_DIR>/build/Linux-x86_64-GCC/softposit.a
         -DAP_INT_INCLUDE_DIR=<PATH_TO_AP_INT_INCLUDE_DIR>
         -DBUILD_UNIT_TEST=1
make
ctest         
```


# Reference 

Marto was first presented at the FPL 2019 conference in Barcelona.
The original article can be found [here](https://hal.archives-ouvertes.fr/hal-02130912v4).

If you use this library, please cite the following reference : 

```Tex
@inproceedings{uguen:hal-02130912,
  TITLE = {{Evaluating the hardware cost of the posit number system}},
  AUTHOR = {Uguen, Yohann and Forget, Luc and de Dinechin, Florent},
  URL = {https://hal.inria.fr/hal-02130912},
  BOOKTITLE = {{FPL 2019 - 29th International Conference on Field-Programmable Logic and Applications (FPL)}},
  ADDRESS = {Barcelona, Spain},
  PAGES = {1-8},
  YEAR = {2019},
  MONTH = Sep,
  PDF = {https://hal.inria.fr/hal-02130912/file/hal_marto_final.pdf},
  HAL_ID = {hal-02130912},
  HAL_VERSION = {v4},
}
```
