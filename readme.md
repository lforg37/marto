To build unit tests:

```
mkdir build 
cd build 
cmake .. -DSOFTPOSIT_H=<PATH_TO_SOFTPOSIT_INSTALL_DIR>/source/include
         -DSOFTPOSIT_LIB=<PATH_TO_SOFTPOSIT_INSTALL_DIR>/build/Linux-x86_64-GCC/softposit.a
         -DAP_INT_INCLUDE_DIR=<PATH_TO_AP_INT_INCLUDE_DIR>
         -DBUILD_UNIT_TEST=1
```