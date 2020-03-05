from Utils import *
from ctypes import *

_Vce0 = 1.8  # V
_Vd0 = 1.3  # V
_Udc = 300  # V
_Toff = 0.32e-6  # sec
_Ton = 0.15e-6  # sec
_Tdead = 3.30e-6  # sec

_Tcomp = 0.0 * (_Tdead + _Ton - _Toff)
# 3.13e-6  #8e-6  # 过补偿 #3.13e-6  # 只补偿死

Rce = 0.04958
Rdiode = 0.05618


class ACIM(object):
    __slot__ = (
        "ial",  # alpha current
        "ibe",  # beta current
        "psi_al",  # stator alpha flux
        "psi_be",  # stator beta flux
        "ual",  # alpha voltage
        "ube",  # beta voltage
        "x0",  #
        "x1",
        "x2",
        "x3",
        "x4",
        "rpm",
        "rpm_cmd",
        "rpm_deriv_cmd",
        "Tload",
        "Tem",
        "Js",
        "npp",
        "mu_m",
        "Ts",
        "Lsigma",
        "rs",
        "rreq",
        "Lmu",
        "Lmu_inv",
        "alpha",
        "Lm",
        "rr",
        "Lls",
        "Llr",
        "Lm_slash_Lr",
        "Lr_slash_Lm",
        "LSigmal",
        "izq",
        "izd",
        "iz",
        "psimq",
        "psimd",
        "psim",
        "im",
        "iqs",
        "ids",
        "iqr",
        "idr",
        "ual_c_dist",
        "ube_c_dist",
        "dist_al",
        "dist_be",
        "RAD_PER_SEC_2_RPM",
        "RPM_2_RAD_PER_SEC",
        "satLUT",
        "satLST",
        "libc",
    )

    def __init__(self):
        # self.ial = 0.0  # Current Alpha Component
        # self.ibe = 0.0  # Current Beta Component
        # self.psi_al = 0.0  #
        # self.psi_be = 0.0  #
        # self.ual = 0.0  # Voltage Alpha Component
        # self.ube = 0.0  # Voltage Beta Component

        # # self.x = [0.0, 0.0, 0.0, 0.0, 0.0]
        # self.x0 = 0.0
        # self.x1 = 0.0
        # self.x2 = 0.0
        # self.x3 = 0.0
        # self.x4 = 0.0  # Motor Status
        # self.rpm = 0.0  # revolutions per minute
        # self.rpm_cmd = 0.0  # speed command
        # self.rpm_deriv_cmd = 0.0  # derivative of speed command
        # self.Tload = 0.0  # Load Torque
        # self.Tem = 0.0  # Electrical Torque

        # self.Js = 0.0  # moment of inertia
        # self.npp = 0.0  # number of pole pairs
        # self.mu_m = 0.0  # =npp/Js
        # self.Ts = 0.0  # sampling time

        # self.Lsigma = 0.0  # Leakage induction
        # self.rs = 0.0  # stator resistance
        # self.rreq = 0.0  # equivalent rotor resistance
        # self.Lmu = 0.0  # magnetic induction
        # self.Lmu_inv = 0.0  # inverse of Lmu
        # self.alpha = 0.0  # inverse of rotor time constant = Lmu/Lr

        # self.Lm = 0.0  # mutual induction
        # self.rr = 0.0  # rotor resistance
        # self.Lls = 0.0  # stator leakage induction
        # self.Llr = 0.0  # rotor leakage induction
        # self.Lm_slash_Lr = 0.0
        # self.Lr_slash_Lm = 0.0

        # self.LSigmal = 0.0

        # self.izq = 0.0
        # self.izd = 0.0
        # self.iz = 0.0

        # self.psimq = 0.0
        # self.psimd = 0.0
        # self.psim = 0.0
        # self.im = 0.0

        # self.iqs = 0.0
        # self.ids = 0.0
        # self.iqr = 0.0
        # self.idr = 0.0

        # self.ual_c_dist = 0.0
        # self.ube_c_dist = 0.0
        # self.dist_al = 0.0
        # self.dist_be = 0.0

        # self.RAD_PER_SEC_2_RPM = 0.0
        # self.RPM_2_RAD_PER_SEC = 0.0
        self.Machine_init()

    def Machine_init(self):
        # self.db = []
        self.x0 = 0.0
        self.x1 = 0.0
        self.x2 = 0.0
        self.x3 = 0.0
        self.x4 = 0.0
        self.rpm = 0.0
        self.rpm_cmd = 0.0
        self.rpm_deriv_cmd = 0.0
        self.Tload = 0.0
        self.Tem = 0.0

        self.Lmu = 0.4482
        # Those parameters are introduced since branch saturation
        self.Lls = 0.0126
        self.Llr = 0.0126
        self.Lm = 0.5 * (
            self.Lmu + np.sqrt(self.Lmu * self.Lmu + 4 * self.Llr * self.Lmu)
        )
        self.Lm_slash_Lr = self.Lm / (self.Lm + self.Llr)
        self.Lr_slash_Lm = (
            self.Lm + self.Llr
        ) / self.Lm  # i_dreq = self.idr*self.Lr_slash_Lm
        self.LSigmal = 1.0 / (1.0 / self.Lls + 1.0 / self.Llr)
        self.Lsigma = self.Lm + self.Lls - self.Lmu  # = Ls * (1.0 - Lm*Lm/Ls/Lr)
        print(
            "Validate: {0} = {1} ?\n".format(
                self.Lsigma,
                (self.Lm + self.Lls)
                * (
                    1.0
                    - self.Lm * self.Lm / (self.Lm + self.Lls) / (self.Lm + self.Llr)
                ),
            )
        )

        self.rreq = 1.69
        self.rs = 3.04
        self.rr = self.rreq * self.Lr_slash_Lm * self.Lr_slash_Lm

        self.alpha = self.rreq / (self.Lmu)
        self.Lmu_inv = 1.0 / self.Lmu

        self.Js = 0.0636  # Awaya92 using im.omg
        self.npp = 2
        self.mu_m = self.npp / self.Js

        self.Ts = MACHINE_TS

        self.ial = 0.0
        self.ibe = 0.0

        self.ual = 0.0
        self.ube = 0.0

        # Those variables are introduced since branch saturation
        self.iqs = 0.0
        self.ids = 0.0
        self.iqr = 0.0
        self.idr = 0.0
        self.psimq = 0.0
        self.psimd = 0.0

        self.satLUT = np.array(pd.read_csv(r"33.txt", header=None))
        self.satLST = (c_float * len(self.satLUT))(*self.satLUT)
        self.libc = cdll.LoadLibrary("./libsatlut.so")
        self.libc.sat_lookup.argtypes = [c_double, POINTER(c_float)]
        self.libc.sat_lookup.restype = c_float

        self.RAD_PER_SEC_2_RPM = 60.0 / (2 * np.pi * self.npp)
        self.RPM_2_RAD_PER_SEC = (2 * np.pi * self.npp) / 60.0

    # def sat_lookup(self, iz, satLUT):
    #     if iz >= (444000) * IZ_STEP:
    #         idx = 444000
    #         return (
    #             (iz - idx * IZ_STEP)
    #             * IZ_STEP_INV
    #             * (satLUT[idx + 1][0] - satLUT[idx][0])
    #             + satLUT[idx][0]
    #         )  # 一阶插值？
    #     # return satLUT[(int)(iz/IZ_STEP+0.5)]  # This yields bad results with obvious oscillation (voltages and currents)
    #     idx = int(iz * IZ_STEP_INV)
    #     return (
    #         (iz - idx * IZ_STEP) * IZ_STEP_INV * (satLUT[idx + 1][0] - satLUT[idx][0])
    #         + satLUT[idx][0]
    #     )  # 一阶插值？

    def collectCurrents(self):
        # Generalised Current by Therrien2013
        self.izq = self.x1 / self.Lls + self.x3 / self.Llr
        self.izd = self.x0 / self.Lls + self.x2 / self.Llr

        self.iz = np.sqrt(self.izd ** 2 + self.izq ** 2)

        if self.iz > 1e-8:
            if SATURATED_MAGNETIC_CIRCUIT:
                self.psim = self.libc.sat_lookup(self.iz, self.satLST)
                self.im = self.iz - self.psim / self.LSigmal
                self.Lm = self.psim / self.im
                self.Lmu = self.Lm * self.Lm / (self.Lm + self.Llr)
                self.Lmu_inv = 1.0 / self.Lmu
                self.alpha = self.rr / (self.Lm + self.Llr)
                self.rreq = self.Lmu * self.alpha
                self.Lsigma = (self.Lls + self.Lm) - self.Lmu
                self.Lm_slash_Lr = self.Lm / (self.Lm + self.Llr)
                self.Lr_slash_Lm = (self.Lm + self.Llr) / self.Lm
            else:
                self.psim = (
                    1.0 / (1.0 / self.Lm + 1.0 / self.Lls + 1.0 / self.Llr) * self.iz
                )
            self.psimq = self.psim / self.iz * self.izq
            self.psimd = self.psim / self.iz * self.izd
        else:
            if selfSIMC_DEBUG:
                print("how to handle zero iz?\n")
            self.psimq = 0.0
            self.psimd = 0.0
        self.iqs = (self.x1 - self.psimq) / self.Lls
        self.ids = (self.x0 - self.psimd) / self.Lls
        self.iqr = (self.x3 - self.psimq) / self.Llr
        self.idr = (self.x2 - self.psimd) / self.Llr
        # # Direct compute is ir from psis psir */
        # self.iqs = (self.x1 - (self.Lm/(self.Lm+self.Llr))*self.x3)/self.Lsigma
        # self.ids = (self.x0 - (self.Lm/(self.Lm+self.Llr))*self.x2)/self.Lsigma
        # self.iqr = (self.x3 - (self.Lm/(self.Lm+self.Lls))*self.x1)/(self.Lm+self.Llr-self.Lm*self.Lm/(self.Lm+self.Lls))
        # self.idr = (self.x2 - (self.Lm/(self.Lm+self.Lls))*self.x0)/(self.Lm+self.Llr-self.Lm*self.Lm/(self.Lm+self.Lls))

    def rK5_satDynamics(self, t, x, fx):
        # argument t is omitted*/

        # STEP ZERO: collect all the currents: is, ir, iz */
        self.collectCurrents()

        # STEP ONE: Inverter Nonlinearity Considered */
        # not CTRL.ual
        # not CTRL.ube

        # STEP TWO: electromagnetic model with full flux states in alpha-beta frame */
        fx[1] = self.ube - self.rs * self.iqs  # 这里是反过来的！Q轴在前面！
        fx[0] = self.ual - self.rs * self.ids
        fx[3] = -self.rr * self.iqr + x[4] * x[2]  # 这里是反过来的！Q轴在前面！
        fx[2] = -self.rr * self.idr - x[4] * x[3]

        # STEP THREE: mechanical model */
        # self.Tem = self.npp*(self.Lm/(self.Lm+self.Llr))*(self.iqs*self.x2-self.ids*self.x3)  # this is not better
        self.Tem = self.npp * (self.iqs * x[0] - self.ids * x[1])
        fx[4] = (self.Tem - self.Tload) * self.mu_m

    def rK555_Sat(self, t, x, hs):
        k1 = [0, 0, 0, 0, 0]
        k2 = [0, 0, 0, 0, 0]
        k3 = [0, 0, 0, 0, 0]
        k4 = [0, 0, 0, 0, 0]
        xk = [0, 0, 0, 0, 0]
        fx = [0, 0, 0, 0, 0]

        self.rK5_satDynamics(t, x, fx)  # timer.t,
        # for i in range(5):
        k1 = [fx[i] * hs for i in range(5)]
        xk = [x[i] + k1[i] * 0.5 for i in range(5)]

        self.rK5_satDynamics(t, xk, fx)  # timer.t+hs/2.,
        # for i in range(5):
        k2 = [fx[i] * hs for i in range(5)]
        xk = [x[i] + k2[i] * 0.5 for i in range(5)]

        self.rK5_satDynamics(t, xk, fx)  # timer.t+hs/2.,
        # for i in range(5):
        k3 = [fx[i] * hs for i in range(5)]
        xk = [x[i] + k3[i] for i in range(5)]

        self.rK5_satDynamics(t, xk, fx)  # timer.t+hs,
        k4 = [fx[i] * hs for i in range(5)]
        tmp = [
            x[i] + (k1[i] + 2 * (k2[i] + k3[i]) + k4[i]) * one_over_six
            for i in range(5)
        ]
        self.x0, self.x1, self.x2, self.x3, self.x4 = tmp
        self.collectCurrents()

    def machine_simulation(self, timebase):
        x = [self.x0, self.x1, self.x2, self.x3, self.x4]
        self.rK555_Sat(timebase, x, self.Ts)
        self.ial = self.ids  # rK555_Sat
        self.ibe = self.iqs  # rK555_Sat
        self.psi_al = self.x2 * self.Lm_slash_Lr  # rK555_Sat
        self.psi_be = self.x3 * self.Lm_slash_Lr  # rK555_Sat

        self.rpm = self.x4 * 60 / (2 * np.pi * self.npp)

        if not np.isnan(self.rpm):
            return False
        else:
            # print("self.rpm is {0}\n", self.rpm)
            return True

    def InverterNonlinearity_SKSul96(self, ual, ube, ial, ibe):
        TM = _Toff - _Ton - _Tdead + _Tcomp  # Sul1996
        Udist = (
            _Udc * TM * TS_INVERSE - _Vce0 - _Vd0
        ) / 6.0  # Udist = (_Udc*TM/1e-4 - _Vce0 - _Vd0) / 6.0
        # Udist = (_Udc*TM*TS_INVERSE) / 6.0
        # Udist = 0.0

        ia = SQRT_2_SLASH_3 * (ial)
        ib = SQRT_2_SLASH_3 * (-0.5 * ial - SIN_DASH_2PI_SLASH_3 * ibe)
        ic = SQRT_2_SLASH_3 * (-0.5 * ial - SIN_2PI_SLASH_3 * ibe)

        # compute in abc frame */
        # ua = SQRT_2_SLASH_3 * (       ual                              )
        # ub = SQRT_2_SLASH_3 * (-0.5 * ual - SIN_DASH_2PI_SLASH_3 * ube )
        # uc = SQRT_2_SLASH_3 * (-0.5 * ual - SIN_2PI_SLASH_3      * ube )
        # ua += Udist * (2*sign(ia) - sign(ib) - sign(ic)) - 0.5*(Rce+Rdiode)*ia
        # ub += Udist * (2*sign(ib) - sign(ic) - sign(ia)) - 0.5*(Rce+Rdiode)*ib
        # uc += Udist * (2*sign(ic) - sign(ia) - sign(ib)) - 0.5*(Rce+Rdiode)*ic
        # UAL_C_DIST = SQRT_2_SLASH_3      * (ua - 0.5*ub - 0.5*uc)  # sqrt(2/3.)
        # UBE_C_DIST = 0.70710678118654746 * (         ub -     uc)  # sqrt(2/3.)*sin(2*pi/3) = sqrt(2/3.)*(sqrt(3)/2)

        # directly compute in alpha-beta frame */
        # CHECK the sign of the distortion voltage!
        # Sul把Udist视为补偿的电压（假定上升下降时间都已经知道了而且是“补偿”上去的）
        self.ual_c_dist = (
            ual
            + np.sqrt(1.5) * Udist * (2 * np.sign(ia) - np.sign(ib) - np.sign(ic))
            - 0.5 * (Rce + Rdiode) * ial
        )
        self.ube_c_dist = (
            ube
            + 3 / np.sqrt(2) * Udist * (np.sign(ib) - np.sign(ic))
            - 0.5 * (Rce + Rdiode) * ibe
        )

    def inverter_model(self, CTRL):

        # 根据给定电压CTRL.ual和实际的电机电流self.ial，计算畸变的逆变器输出电压self.ual。
        if INVERTER_NONLINEARITY:

            self.InverterNonlinearity_SKSul96(CTRL.ual, CTRL.ube, self.ial, self.ibe)
            # InverterNonlinearity_Tsuji01
            self.ual = self.ual_c_dist
            self.ube = self.ube_c_dist

            # Distorted voltage (for visualization purpose)
            self.dist_al = self.ual - CTRL.ual
            self.dist_be = self.ube - CTRL.ube
        else:
            self.ual = CTRL.ual
            self.ube = CTRL.ube
