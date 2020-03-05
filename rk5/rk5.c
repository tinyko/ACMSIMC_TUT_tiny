#include "rk5.h"

struct args_in ai;
struct args_out ao;

int collectCurrents()
{
    double im, psimq, psimd, alpha;
    ai.izq = ai.x1 / ai.Lls + ai.x3 / ai.Llr;
    ai.izd = ai.x0 / ai.Lls + ai.x2 / ai.Llr;

    ai.iz = sqrt(ai.izd * ai.izd + ai.izq * ai.izq);

    if (ai.iz > 1e-8)
    {
        if (SATURATED_MAGNETIC_CIRCUIT)
        {
            ai.psim = sat_lookup(ai.iz, satLUT);
            im = ai.iz - ai.psim / ai.LSigmal;
            ai.Lm = ai.psim / im;
            ai.Lmu = ai.Lm * ai.Lm / (ai.Lm + ai.Llr);
            ai.Lmu_inv = 1.0 / ai.Lmu;
            ai.alpha = ai.rr / (ai.Lm + ai.Llr);
            ai.rreq = ai.Lmu * ai.alpha;
            ai.Lsigma = (ai.Lls + ai.Lm) - ai.Lmu;
            ai.Lm_slash_Lr = ai.Lm / (ai.Lm + ai.Llr);
            ai.Lr_slash_Lm = (ai.Lm + ai.Llr) / ai.Lm;
        }
        else
        {
            ai.psim = (1.0 / (1.0 / ai.Lm + 1.0 / ai.Lls + 1.0 / ai.Llr) * ai.iz);
        }
        psimq = ai.psim / ai.iz * ai.izq;
        psimd = ai.psim / ai.iz * ai.izd;
    }
    else
    {
        psimq = 0.0;
        psimd = 0.0;
    }

    ai.iqs = (ai.x1 - psimq) / ai.Lls;
    ai.ids = (ai.x0 - psimd) / ai.Lls;
    ai.iqr = (ai.x3 - psimq) / ai.Llr;
    ai.idr = (ai.x2 - psimd) / ai.Llr;
}

int rK5_satDynamics(double t, double x[], double fx[])
{
    collectCurrents();

    // # STEP ONE: Inverter Nonlinearity Considered */
    // # not CTRL.ual
    // # not CTRL.ube

    // # STEP TWO: electromagnetic model with full flux states in alpha-beta frame */
    fx[1] = ai.ube - ai.rs * ai.iqs; //# 这里是反过来的！Q轴在前面！
    fx[0] = ai.ual - ai.rs * ai.ids;
    fx[3] = -ai.rr * ai.iqr + x[4] * x[2]; //# 这里是反过来的！Q轴在前面！
    fx[2] = -ai.rr * ai.idr - x[4] * x[3];

    // # STEP THREE: mechanical model */
    // # ai.Tem = ai.npp*(ai.Lm/(ai.Lm+ai.Llr))*(ai.iqs*ai.x2-ai.ids*ai.x3)  # this is not better
    ai.Tem = ai.npp * (ai.iqs * x[0] - ai.ids * x[1]);
    fx[4] = (ai.Tem - ai.Tload) * ai.mu_m;
}

int copy_args_out()
{
    ao.x0 = ai.x0;
    ao.x1 = ai.x1;
    ao.x2 = ai.x2;
    ao.x3 = ai.x3;
    ao.x4 = ai.x4;                   //IO
    ao.iqs = ai.iqs;                 //IO
    ao.ids = ai.ids;                 //IO
    ao.iqr = ai.iqr;                 //IO
    ao.idr = ai.idr;                 //IO
    ao.Tem = ai.Tem;                 //IO
    ao.izq = ai.izq;                 //IO
    ao.izd = ai.izd;                 //IO
    ao.iz = ai.iz;                   //IO
    ao.psim = ai.psim;               //IO
    ao.Lm = ai.Lm;                   //IO
    ao.LSigmal = ai.LSigmal;         //IO
    ao.Lmu = ai.Lmu;                 //IO
    ao.Lmu_inv = ai.Lmu_inv;         //IO
    ao.alpha = ai.alpha;             //IO
    ao.rreq = ai.rreq;               //IO
    ao.Lsigma = ai.Lsigma;           //IO
    ao.Lr_slash_Lm = ai.Lr_slash_Lm; //IO
    ao.Lm_slash_Lr = ai.Lm_slash_Lr; //IO
}

struct args_out rK555_Sat(struct args_in ai_in)
{
    ai = ai_in;
    int i;
    double t = ai.timebase;
    double hs = ai.Ts;
    double k1[5], k2[5], k3[5], k4[5], xk[5], fx[5] = {0};
    double x[] = {ai.x0, ai.x1, ai.x2, ai.x3, ai.x4};

    rK5_satDynamics(t, x, fx); //# timer.t,
    for (i = 0; i < 5; i++)
    {
        k1[i] = fx[i] * hs;
        xk[i] = x[i] + k1[i] * 0.5;
    }

    rK5_satDynamics(t, xk, fx); //# timer.t+hs/2.,
    for (i = 0; i < 5; i++)
    {
        k2[i] = fx[i] * hs;
        xk[i] = x[i] + k2[i] * 0.5;
    }

    rK5_satDynamics(t, xk, fx); //# timer.t+hs/2.,
    for (i = 0; i < 5; i++)
    {
        k3[i] = fx[i] * hs;
        xk[i] = x[i] + k3[i];
    }

    rK5_satDynamics(t, xk, fx); //# timer.t+hs,
    for (i = 0; i < 5; i++)
    {
        k4[i] = fx[i] * hs;
        x[i] = x[i] + (k1[i] + 2 * (k2[i] + k3[i]) + k4[i]) / 6;
    }
    ai.x0 = x[0];
    ai.x1 = x[1];
    ai.x2 = x[2];
    ai.x3 = x[3];
    ai.x4 = x[4];
    collectCurrents();
    copy_args_out();
    return ao;
}
