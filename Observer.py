from Utils import MACHINE_TYPE, INDUCTION_MACHINE, TS, np


class Tajima(object):
    def __init__(self):
        self.K_PEM = 0.0
        self.omega_syn = 0.0
        self.e_M = 0.0
        self.omega_sl = 0.0
        self.omg = 0.0


class InductionMachine(object):
    def __init__(self):
        self.zeros = np.zeros(2)
        self.u_s = self.zeros
        self.i_s = self.zeros
        self.us_curr = self.zeros
        self.is_curr = self.zeros
        self.us_prev = self.zeros
        self.is_prev = self.zeros

        self.Js = 0.0
        self.Js_inv = 0.0

        self.rs = 0.0
        self.rreq = 0.0

        self.alpha = 0.0
        self.Lsigma = 0.0
        self.Lsigma_inv = 0.0
        self.Lmu = 0.0
        self.Lmu_inv = 0.0

        self.npp = 0.0
        self.omg = 0.0
        self.theta_r = 0.0
        self.theta_d = 0.0


class Observer(object):
    def __init__(self, IM):
        self.tajima = Tajima()
        self.im = InductionMachine()
        self.acm_init(IM)
        self.ob_init()

    def acm_init(self, IM):
        self.im.u_s = self.im.zeros
        self.im.i_s = self.im.zeros
        self.im.us_curr = self.im.zeros
        self.im.is_curr = self.im.zeros
        self.im.us_prev = self.im.zeros
        self.im.is_prev = self.im.zeros

        self.im.Js = IM.Js
        self.im.Js_inv = 1.0 / self.im.Js
        self.im.rs = IM.rs
        self.im.rreq = IM.rreq
        self.im.alpha = IM.alpha
        self.im.Lsigma = IM.Lsigma
        self.im.Lsigma_inv = 1 / self.im.Lsigma
        self.im.Lmu = IM.Lmu
        self.im.Lmu_inv = IM.Lmu_inv
        self.im.npp = IM.npp
        self.im.omg = 0
        self.im.theta_r = 0.0

    def ob_init(self):
        self.Js = self.im.Js
        self.Js_inv = self.im.Js_inv
        self.rs = self.im.rs
        self.rreq = self.im.rreq
        self.alpha = self.im.alpha
        self.Lsigma = self.im.Lsigma
        self.Lsigma_inv = self.im.Lsigma_inv
        self.Lmu = self.im.Lmu
        self.Lmu_inv = self.im.Lmu_inv
        self.npp = self.im.npp
        self.omg = self.im.omg
        self.im.theta_r = 0.0
        self.im.theta_d = 0.0

        self.psi_mu_al = 0.0
        self.psi_mu_be = 0.0
        self.psi_s_al = 0.0
        self.psi_s_be = 0.0

    def measurement(self, IM, CTRL):
        self.im.us_curr[0] = CTRL.ual
        self.im.us_curr[1] = CTRL.ube
        self.im.us_prev[0] = self.im.us_curr[0]
        self.im.us_prev[1] = self.im.us_curr[1]

        if MACHINE_TYPE == INDUCTION_MACHINE:
            self.im.is_curr[0] = IM.ial
            self.im.is_curr[1] = IM.ibe
            self.im.omg = IM.x4
            # self.im.theta_r = IM.x5

    def observation(self, CTRL):
        rotor_flux_cmd = CTRL.rotor_flux_cmd
        iMs = CTRL.iMs
        iTs = CTRL.iTs
        uMs_cmd = CTRL.uMs_cmd
        uTs_cmd = CTRL.uTs_cmd

        # Speed estimation: Tajima1996
        self.tajima.K_PEM = 0.1  # 0.1 ~ 2
        self.tajima.omega_syn = (uTs_cmd - self.rs * iTs) / (
            rotor_flux_cmd + self.Lsigma * iMs
        )
        self.tajima.e_M = (
            uMs_cmd - self.rs * iMs + self.tajima.omega_syn * self.Lsigma * iTs
        )
        self.tajima.omega_syn -= self.tajima.K_PEM * self.tajima.e_M

        self.tajima.omega_sl = self.rreq * iTs / rotor_flux_cmd
        self.tajima.omg = (
            self.tajima.omega_syn - self.tajima.omega_sl
        )  # Instantaneous Velocity Computation

        # Flux estimation 1: Voltage model
        # (implemented by shitty Euler method for now)

        deriv_psi_s_al = self.im.us_curr[0] - self.rs * self.im.is_curr[0]
        deriv_psi_s_be = self.im.us_curr[1] - self.rs * self.im.is_curr[1]
        # psi_s_al = psi_s_be = 0

        self.psi_s_al += TS * deriv_psi_s_al
        self.psi_s_be += TS * deriv_psi_s_be
        self.psi_mu_al = self.psi_s_al - self.Lsigma * self.im.i_s[0]
        self.psi_mu_be = self.psi_s_be - self.Lsigma * self.im.i_s[1]

        # Flux estimation 2: Current model (this is a bad flux estimator)
