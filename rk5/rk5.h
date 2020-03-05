#include <stdio.h>
#include <math.h>
#include "satlut.h"
#ifndef SATURATED_MAGNETIC_CIRCUIT
#define SATURATED_MAGNETIC_CIRCUIT 1
#endif
#define one_over_six (double)(1 / 6)
struct args_in
{
    double Llr;      //I
    double Lls;      //I
    double timebase; //I
    double Ts;       //I
    double ube;      //I
    double ual;      //I
    double rs;       //I
    double rr;       //I
    double npp;      //I
    double Tload;    //I
    double mu_m;     //I
    double x0;
    double x1;
    double x2;
    double x3;
    double x4;          //IO
    double iqs;         //IO
    double ids;         //IO
    double iqr;         //IO
    double idr;         //IO
    double Tem;         //IO
    double izq;         //IO
    double izd;         //IO
    double iz;          //IO
    double psim;        //IO
    double Lm;          //IO
    double LSigmal;     //IO
    double Lmu;         //IO
    double Lmu_inv;     //IO
    double alpha;       //IO
    double rreq;        //IO
    double Lsigma;      //IO
    double Lm_slash_Lr; //IO
    double Lr_slash_Lm; //IO
};

struct args_out
{
    double x0;
    double x1;
    double x2;
    double x3;
    double x4;          //IO
    double iqs;         //IO
    double ids;         //IO
    double iqr;         //IO
    double idr;         //IO
    double Tem;         //IO
    double izq;         //IO
    double izd;         //IO
    double iz;          //IO
    double psim;        //IO
    double Lm;          //IO
    double LSigmal;     //IO
    double Lmu;         //IO
    double Lmu_inv;     //IO
    double alpha;       //IO
    double rreq;        //IO
    double Lsigma;      //IO
    double Lm_slash_Lr; //IO
    double Lr_slash_Lm; //IO
};

int collectCurrent();
int rK5_satDynamics(double t, double x[], double fx[]);
struct args_out rK555_Sat(struct args_in ai);
int copy_args_out();