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


class args_out(Structure):
    _fields_ = [
        ("x0", c_double),
        ("x1", c_double),
        ("x2", c_double),
        ("x3", c_double),
        ("x4", c_double),
        ("iqs", c_double),
        ("ids", c_double),
        ("iqr", c_double),
        ("idr", c_double),
        ("Tem", c_double),
        ("izq", c_double),
        ("izd", c_double),
        ("iz", c_double),
        ("psim", c_double),
        ("Lm", c_double),
        ("LSigmal", c_double),
        ("Lmu", c_double),
        ("Lmu_inv", c_double),
        ("alpha", c_double),
        ("rreq", c_double),
        ("Lsigma", c_double),
        ("Lm_slash_Lr", c_double),
        ("Lr_slash_Lm", c_double),
    ]


class args_in(Structure):
    _fields_ = [
        ("Llr", c_double),
        ("Lls", c_double),
        ("timebase", c_double),
        ("Ts", c_double),
        ("ube", c_double),
        ("ual", c_double),
        ("rs", c_double),
        ("rr", c_double),
        ("npp", c_double),
        ("Tload", c_double),
        ("mu_m", c_double),
        ("x0", c_double),
        ("x1", c_double),
        ("x2", c_double),
        ("x3", c_double),
        ("x4", c_double),
        ("iqs", c_double),
        ("ids", c_double),
        ("iqr", c_double),
        ("idr", c_double),
        ("Tem", c_double),
        ("izq", c_double),
        ("izd", c_double),
        ("iz", c_double),
        ("psim", c_double),
        ("Lm", c_double),
        ("LSigmal", c_double),
        ("Lmu", c_double),
        ("Lmu_inv", c_double),
        ("alpha", c_double),
        ("rreq", c_double),
        ("Lsigma", c_double),
        ("Lm_slash_Lr", c_double),
        ("Lr_slash_Lm", c_double),
    ]


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

        self.izq = 0.0
        self.izd = 0.0
        self.iz = 0.0
        self.psim = 0.0

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

        self.ao_out = args_out()
        self.librk4 = cdll.LoadLibrary("./rk4/librk4.so")
        self.librk4.rK4_Sat.argtypes = [args_in]
        self.librk4.rK4_Sat.restype = args_out
        self.ai_in = args_in()
        self.ao_out = args_out()
        self.RAD_PER_SEC_2_RPM = 60.0 / (2 * np.pi * self.npp)
        self.RPM_2_RAD_PER_SEC = (2 * np.pi * self.npp) / 60.0

    def copy_vars_to_in(self):
        self.ai_in.Llr = self.Llr
        self.ai_in.Lls = self.Lls
        self.ai_in.timebase = self.timebase
        self.ai_in.Ts = self.Ts
        self.ai_in.ube = self.ube
        self.ai_in.ual = self.ual
        self.ai_in.rs = self.rs
        self.ai_in.rr = self.rr
        self.ai_in.npp = self.npp
        self.ai_in.Tload = self.Tload
        self.ai_in.mu_m = self.mu_m
        self.ai_in.x0 = self.x0
        self.ai_in.x1 = self.x1
        self.ai_in.x2 = self.x2
        self.ai_in.x3 = self.x3
        self.ai_in.x4 = self.x4
        self.ai_in.iqs = self.iqs
        self.ai_in.ids = self.ids
        self.ai_in.iqr = self.iqr
        self.ai_in.idr = self.idr
        self.ai_in.Tem = self.Tem
        self.ai_in.izq = self.izq
        self.ai_in.izd = self.izd
        self.ai_in.iz = self.iz
        self.ai_in.psim = self.psim
        self.ai_in.Lm = self.Lm
        self.ai_in.LSigmal = self.LSigmal
        self.ai_in.Lmu = self.Lmu
        self.ai_in.Lmu_inv = self.Lmu_inv
        self.ai_in.alpha = self.alpha
        self.ai_in.rreq = self.rreq
        self.ai_in.Lsigma = self.Lsigma
        self.ai_in.Lr_slash_Lm = self.Lr_slash_Lm
        self.ai_in.Lm_slash_Lr = self.Lm_slash_Lr

    def copy_vars_from_out(self):
        self.x0 = self.ao_out.x0
        self.x1 = self.ao_out.x1
        self.x2 = self.ao_out.x2
        self.x3 = self.ao_out.x3
        self.x4 = self.ao_out.x4
        self.iqs = self.ao_out.iqs
        self.ids = self.ao_out.ids
        self.iqr = self.ao_out.iqr
        self.idr = self.ao_out.idr
        self.Tem = self.ao_out.Tem
        self.izq = self.ao_out.izq
        self.izd = self.ao_out.izd
        self.iz = self.ao_out.iz
        self.psim = self.ao_out.psim
        self.Lm = self.ao_out.Lm
        self.LSigmal = self.ao_out.LSigmal
        self.Lmu = self.ao_out.Lmu
        self.Lmu_inv = self.ao_out.Lmu_inv
        self.alpha = self.ao_out.alpha
        self.rreq = self.ao_out.rreq
        self.Lsigma = self.ao_out.Lsigma
        self.Lr_slash_Lm = self.ao_out.Lr_slash_Lm
        self.Lm_slash_Lr = self.ao_out.Lm_slash_Lr

    def machine_simulation(self, timebase):
        self.timebase = timebase
        self.copy_vars_to_in()
        self.ao_out = self.librk4.rK4_Sat(self.ai_in)
        self.copy_vars_from_out()
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
