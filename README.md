# Tiny Digital Filter Library readme


Yorck Ramachers (Warwick)

Last updated March 4, 2023

## Description

Depends on: ROOT, tested with version 6.14.xx

To build it, do

``` console
$ ls
CMakeLists.txt  DFilter.cpp  DFilter.h  LICENSE  README.md  test

$ mkdir build
$ cd build
$ cmake ..
...
$ make
...
... If you are developing the module, you can test it by calling
$ make test
...
... or obtain more detail on the tests and launch in the build directory
$ ctest -V
$ ./testdf
```

The build will create the `DFilter.so` shared library and (currently)
two test executables in the build directory (from which you can run
the executables, no problem).

The last testing command is a test executable outside
ctest since it creates a ROOT output file for checking
rather than require checks from the catch library.

## Content:
- A Low-pass Butterworth filter object of any order.
- A matched filter object, specialized on finding a small template in a much longer data set.
- A moving average object - smoothing while preserving sharp features.
- A low-pass object - classic RC filter
- A high-pass object - classic CR filter

Check the testdfilter.cpp source in the test directory on usage.
