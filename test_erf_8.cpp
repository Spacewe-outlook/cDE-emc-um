#include <iostream>
#include <fstream>
#include <math.h>
#include<cstdlib>
#include<ctime>
using namespace std;
#define N  99999 
#define PI 3.1415926535897932384626433832795029

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
double myerfinv03(double x){
    double temp,re;
    temp = x+0.261799388*x*x*x+0.143931731*x*x*x*x*x+0.097663620*x*x*x*x*x*x*x+0.073299079*x*x*x*x*x*x*x*x*x+0.058372501*x*x*x*x*x*x*x*x*x*x*x+0.048336063*x*x*x*x*x*x*x*x*x*x*x*x*x+0.041147395*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x;
    re = temp*(sqrt(PI)/2);
    return re;

}
double myerf01(double x){
    double a1=0.278393,a2=0.230389,a3=0.000972,a4=0.078108;
    double t,rer;
    double txx,ttt;
    txx=x*x;
    t= 1+a1*x+a2*txx+a3*txx*x+a4*txx*txx;
    ttt = t*t;
    rer=1-1/(ttt*ttt);

    return rer;
}
double myerf02(double x){
    double a1=0.3480242,a2=-0.0958798,a3=0.7478556;
    double t,rer,t22;
    t=1/(1+0.47047*x);
    t22 =t*t;
    rer = 1-(a1*t+a2*t22+a3*t22*t)*exp(-x*x);
    if(x<0){
        rer = -rer;
    }
    return rer;
}
double myerf03(double x){
    double a1=0.0705230784,a2=0.0422820123,a3=0.0092705272,a4= 0.0001520143,a5=0.0002765672,a6=0.0000430638;
    double t,rer,t33,t88,t44,t22;
    t33=x*x*x;
    t=1+a1*x+a2*x*x+a3*t33+a4*t33*x+a5*t33*x*x+a6*t33*t33;
    t22=t*t;
    t44=t22*t22;
    t88=t44*t44;
    rer = 1-1/(t88*t88);
    if(x<0){
        rer = -rer;
    }
    return rer;
}
double myerf04(double x){
    double a1= 0.254829592,a2=-0.284496736,a3=1.421413741,a4=-1.453152027,a5=1.061405429;
    double t,rer,t33,t22;
    t=1/(1+0.3275911*x);
    t22 = t*t;
    t33 = t22*t;
    rer = 1-(a1*t+a2*t22+a3*t33+a4*t22*t22+a5*t33*t22)*exp(-x*x);
    if(x<0){
        rer = -rer;
    }
    return rer;
}


double generateCDFInv_0(double r,double mu,double sigma){
    double erfA,erfB,erftemp,samplerand;
    erftemp = sqrt(2)*sigma;
    erfA = erf((mu+1)/erftemp);
    erfB = erf((mu-1)/erftemp);
    samplerand = myerfinv(-erfA-r*erfB+r*erfA)*erftemp+mu;
    return samplerand;
}
double generateCDFInv_01(double r,double mu,double sigma){
    double erfA,erfB,erftemp,samplerand;
    erftemp = sqrt(2)*sigma;
    erfA = myerf01((mu+1)/erftemp);
    erfB = myerf01((mu-1)/erftemp);
    samplerand = myerfinv(-erfA-r*erfB+r*erfA)*erftemp+mu;
    return samplerand;
}
double generateCDFInv_02(double r,double mu,double sigma){
    double erfA,erfB,erftemp,samplerand;
    erftemp = sqrt(2)*sigma;
    erfA = myerf02((mu+1)/erftemp);
    erfB = myerf02((mu-1)/erftemp);
    samplerand = myerfinv(-erfA-r*erfB+r*erfA)*erftemp+mu;
    return samplerand;
}
double generateCDFInv_03(double r,double mu,double sigma){
    double erfA,erfB,erftemp,samplerand;
    erftemp = sqrt(2)*sigma;
    erfA = myerf03((mu+1)/erftemp);
    erfB = myerf03((mu-1)/erftemp);
    samplerand = myerfinv(-erfA-r*erfB+r*erfA)*erftemp+mu;
    return samplerand;
}
double generateCDFInv_04(double r,double mu,double sigma){
    double erfA,erfB,erftemp,samplerand;
    erftemp = sqrt(2)*sigma;
    erfA = myerf04((mu+1)/erftemp);
    erfB = myerf04((mu-1)/erftemp);
    samplerand = myerfinv(-erfA-r*erfB+r*erfA)*erftemp+mu;
    return samplerand;
}

double generateCDFInv_11(double r,double mu,double sigma){
    double erfA,erfB,erftemp,samplerand;
    erftemp = sqrt(2)*sigma;
    erfA = myerf01((mu+1)/erftemp);
    erfB = myerf01((mu-1)/erftemp);
    samplerand = myerfinv03(-erfA-r*erfB+r*erfA)*erftemp+mu;
    return samplerand;
}
double generateCDFInv_12(double r,double mu,double sigma){
    double erfA,erfB,erftemp,samplerand;
    erftemp = sqrt(2)*sigma;
    erfA = myerf02((mu+1)/erftemp);
    erfB = myerf02((mu-1)/erftemp);
    samplerand = myerfinv03(-erfA-r*erfB+r*erfA)*erftemp+mu;
    return samplerand;
}
double generateCDFInv_13(double r,double mu,double sigma){
    double erfA,erfB,erftemp,samplerand;
    erftemp = sqrt(2)*sigma;
    erfA = myerf03((mu+1)/erftemp);
    erfB = myerf03((mu-1)/erftemp);
    samplerand = myerfinv03(-erfA-r*erfB+r*erfA)*erftemp+mu;
    return samplerand;
}
double generateCDFInv_14(double r,double mu,double sigma){
    double erfA,erfB,erftemp,samplerand;
    erftemp = sqrt(2)*sigma;
    erfA = myerf04((mu+1)/erftemp);
    erfB = myerf04((mu-1)/erftemp);
    samplerand = myerfinv03(-erfA-r*erfB+r*erfA)*erftemp+mu;
    return samplerand;
}
// test erf
// int main(){
    
//     // int Func_all_num = 8;
//     // for(int i=0;i<Func_all_num;i++){
//         int GenNUM = 1000;
//         printf("begin %d\n",GenNUM);
//         srand(time(NULL));
//         double err_01 = 0;
//         double err_02 = 0;
//         double err_03 = 0;
//         double err_04 = 0;
//         double  or_A=0,or_1=0,or_2=0,or_3=0,or_4=0,or_5=0,or_6=0,or_7=0,or_8=0,rand_r=0;
//         for(int i=0;i<GenNUM;i++){
//             rand_r = rand() % (N + 1) / (float)(N + 1);
//             or_A = erf(rand_r);
//             or_1 = myerf01(rand_r);
//             or_2 = myerf01(rand_r);
//             or_3 = myerf01(rand_r);
//             or_4 = myerf01(rand_r);
//             err_01 = err_01 + abs(or_A-or_1);
//             err_02 = err_02 + abs(or_A-or_2);
//             err_03 = err_03 + abs(or_A-or_3);
//             err_04 = err_04 + abs(or_A-or_4);
//             cout << rand_r <<","<<err_01<<","<<err_02<<","<<err_03<<","<<err_04<<endl; 


//         }
//         err_01 = err_01/GenNUM;
//         err_02 = err_02/GenNUM;
//         err_03 = err_03/GenNUM;
//         err_04 = err_04/GenNUM;
//         cout <<err_01<<","<<err_02<<","<<err_03<<","<<err_04<<endl; 


// }

int main(){
    
        int GenNUM = 1000000;
        printf("begin %d\n",GenNUM);
        srand(time(NULL));
        double err_01 = 0;
        double err_02 = 0;
        double err_03 = 0;
        double err_04 = 0;
        double err_11 = 0;
        double err_12 = 0;
        double err_13 = 0;
        double err_14 = 0;
        double  or_A=0,or_1=0,or_2=0,or_3=0,or_4=0,or_5=0,or_6=0,or_7=0,or_8=0,rand_r=0;
        for(int i=0;i<GenNUM;i++){
            rand_r = rand() % (N + 1) / (float)(N + 1);
            or_A = generateCDFInv_0(rand_r,0,10);
            or_1 = generateCDFInv_01(rand_r,0,10);
            or_2 = generateCDFInv_02(rand_r,0,10);
            or_3 = generateCDFInv_03(rand_r,0,10);
            or_4 = generateCDFInv_04(rand_r,0,10);
            or_5 = generateCDFInv_11(rand_r,0,10);
            or_6 = generateCDFInv_12(rand_r,0,10);
            or_7 = generateCDFInv_13(rand_r,0,10);
            or_8 = generateCDFInv_14(rand_r,0,10);
            err_01 = err_01 + abs(or_A-or_1);
            err_02 = err_02 + abs(or_A-or_2);
            err_03 = err_03 + abs(or_A-or_3);
            err_04 = err_04 + abs(or_A-or_4);
            err_11 = err_11 + abs(or_A-or_5);
            err_12 = err_12 + abs(or_A-or_6);
            err_13 = err_13 + abs(or_A-or_7);
            err_14 = err_14 + abs(or_A-or_8);
            cout << rand_r <<","<<err_01<<","<<err_02<<","<<err_03<<","<<err_04<<endl; 


        }
        err_01 = err_01/GenNUM;
        err_02 = err_02/GenNUM;
        err_03 = err_03/GenNUM;
        err_04 = err_04/GenNUM;
        err_11 = err_11/GenNUM;
        err_12 = err_12/GenNUM;
        err_13 = err_13/GenNUM;
        err_14 = err_14/GenNUM;
        cout <<err_01<<","<<err_02<<","<<err_03<<","<<err_04<<","<<err_11<<","<<err_12<<","<<err_13<<","<<err_14<<endl; 


}