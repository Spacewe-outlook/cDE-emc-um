# cDE-emc-um

> The decompression password of the .zip file is the paper ID "TII-****". 
> After the paper is published, it is completely open source, and the decompression password is released.

1. The cDE-emc-um algorithm is contained in the `cDE-emc-um.zip` file and can be run directly using C++. 

2. The `cDE-emc-um for two benchmatk function in FPGA.zip` file contains the FPGA implementation of the cDE-emc-um algorithm for different benchmark test functions, including the hardware implementation of the cDE-emc-um algorithm. The `.sdk` folder contains the program in ARM that calls cDE-emc-um, UnSS, and UnOF. 

3. The `cDE-emc-um for two problems in FPGA.zip` file contains a unified implementation of the cDE-emc-um algorithm for two different problems, including the hardware implementation of UnSS and UnOF. The `.sdk` file implements the unified packaging of data, which processes different problems in a unified manner. 

## Implementation of $erf$ and $erf^{-1}$
Ref : Abramowitz, M. and Stegun, I. A. (Eds.). "Error Function and Fresnel Integrals." Ch. 7 in Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 297-309, 1972.

The mathematical representation of $erf$ and $erf^{-1}$ is shown below.

<img width="266" alt="image" src="https://github.com/Spacewe-outlook/cDE-emc-um/assets/63773526/f31fc33c-5681-4361-88fa-8f78455ffff3">

The reference site is https://mathworld.wolfram.com/Erf.html.

<img width="297" alt="image" src="https://github.com/Spacewe-outlook/cDE-emc-um/assets/63773526/44200e57-65c7-4466-88ee-b654b963c859">

The reference site is https://mathworld.wolfram.com/InverseErf.html. 

The generic $erf$  function can be called from the C++ implementation of the  $erf$ function by using the following C++ header file

```
#include <ap_int.h>
#include <string.h>
#include <ap_fixed.h>
#include "ap_axi_sdata.h"
```

But $erf^{-1}$ is not implemented in C++, so it needs to be reproduced. The generic $erf^{-1}$ function is implemented with the following C++ code

```
double myerfinv(double yo){
    double z,xo = 0;
    double acoef[4] = { 0.886226899,-1.645349621, 0.914624893, -0.140543331};
    double bcoef[4] = {-2.118377725,  1.442710462, -0.329097515,  0.012229801};
    double ccoef[4] = {-1.970840454, -1.624906493,  3.429567803,  1.641345311};
    double dcoef[2] = { 3.543889200,  1.637067800};
    double y0 = 0.7;
    if (abs(yo) <= y0){
        z = yo*yo;
        xo = yo* (((acoef[3]*z+acoef[2])*z+acoef[1])*z+acoef[0])/((((bcoef[3]*z+bcoef[2])*z+bcoef[1])*z+bcoef[0])*z+1);
    }
    if(( y0 < yo ) && (yo <  1)){
        z = sqrt(-log((1-yo)/2));
        xo = (((ccoef[3]*z+ccoef[2])*z+ccoef[1])*z+ccoef[0])/((dcoef[1]*z+dcoef[0])*z+1);
    }
    if((-y0 > yo ) && (yo > -1)){
        z = sqrt(-log((1+yo)/2));
        xo = -(((ccoef[3]*z+ccoef[2])*z+ccoef[1])*z+ccoef[0]) / ((dcoef[1]*z+dcoef[0])*z+1);
    }
    return xo;
}
```
