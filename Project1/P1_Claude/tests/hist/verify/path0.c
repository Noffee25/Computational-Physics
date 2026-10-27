#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int argc,char**argv){ for(int a=1;a<argc;a++){ long long seed=atoll(argv[a]);
 double y=fmod((double)seed*0.6180339887498949,1.0), p2=-1,p1=-1; long long n=0;
 while(y!=0.0 && n<200000000LL){ p2=p1; p1=y; y=4.0*y*(1.0-y); n++; }
 printf("%lld steps %lld prev2 %.17g prev1 %.17g now %.17g\n",seed,n,p2,p1,y);} return 0;}
