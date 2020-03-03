from Utils import *
import Macros as mc


class s_curve(object):
    __slot__ = (
        "V_T",
        "ACC",
        "DEC",
        "Start",
        "dt",
        "Slow",
        "F1",
        "F2",
        "F3",
        "F4",
        "F5",
        "Target",
        "Stat",
        "V_O",
        "LowV",
        "Busy",
        "acc",
        "dec",
    )

    def __init__(self):
        self.init()

    def init(self):
        self.V_T = 0.0
        self.ACC = 2
        self.DEC = 2
        self.Start = True
        self.dt = 0.125  # micro second
        self.Slow = False
        self.F1 = 0.0
        self.F2 = 0.0
        self.F3 = 0.0
        self.F4 = 0.0
        self.F5 = 0.0
        self.Target = 0.0
        self.Stat = 0
        self.V_O = 0.0
        self.LowV = 2
        self.Busy = False
        self.acc = False
        self.dec = True

    def ref(self):
        if not self.Busy:
            if self.Slow and self.Start:
                self.Target = self.LowV
                self.Stat = 1
            elif not self.Start:
                self.Target = 0
                self.Stat = 2
            else:
                self.Target = self.V_T
                self.Stat = 0

        if (np.abs(self.V_O - self.Target) > 0.001) and (not self.Busy):
            if self.Target > self.V_O:
                h = self.dt / ((self.Target - self.V_O) / 60 / self.ACC) / 1000
                acc = 1
            elif self.Target < self.V_O:
                h = self.dt / ((self.V_O - self.Target) / 60 / self.DEC) / 1000
                dec = 1
            h3 = h * h * h
            h4 = h3 * h
            h5 = h4 * h
            Initial = self.V_O
            A = -6 * Initial + 6 * self.Target
            B = 15 * Initial - 15 * self.Target
            C = -10 * Initial + 10 * self.Target
            self.F5 = (121 * A * h5) / 16 + 5 * B * h4 + (13 * C * h3) / 4
            self.F4 = (165 * A * h5) / 2 + 29 * B * h4 + 9 * C * h3
            self.F3 = 255 * A * h5 + 48 * B * h4 + 6 * C * h3
            self.F2 = 300 * A * h5 + 24 * B * h4
            self.F1 = 120 * A * h5
            self.Busy = True
        elif (np.abs(self.V_O - self.Target) > 0.001) and self.Busy:
            self.V_O = self.V_O + self.F5
            self.F5 = self.F5 + self.F4
            self.F4 = self.F4 + self.F3
            self.F3 = self.F3 + self.F2
            self.F2 = self.F2 + self.F1
            if (self.acc) and (self.V_O >= self.Target):
                self.V_O = self.Target
                self.Busy = False
                self.acc = False
            elif (self.dec) and (self.V_O <= self.Target):
                self.V_O = self.Target
                self.Busy = False
                self.dec = False
        else:
            self.V_O = self.Target
            self.Busy = False
            self.acc = False
            self.dec = False

    def speed_ref(self, timebase, IM):
        if timebase > 5.0 and timebase <= 10.0:
            self.V_T = 100
            self.ref()
            tmp = self.V_O
        elif timebase > 10.0:
            self.V_T = -100
            self.ref()
            tmp = self.V_O
        else:
            self.ref()
            tmp = self.V_O
        IM.rpm_cmd = tmp


class PI_REG:
    __slot__ = ("Kp", "Ti", "Ki", "i_state", "i_limit")

    def __init__(self):
        self.Kp = 0.0
        self.Ti = 0.0
        self.Ki = 0.0  # Ki = Kp/Ti*TS
        self.i_state = 0.0
        self.i_limit = 0.0


class CTRL0(object):
    __slot__ = (
        "vc_count",
        "timebase",
        "ual",
        "ube",
        "rs",
        "rreq",
        "Lsigma",
        "alpha",
        "Lmu",
        "Lmu_inv",
        "Tload",
        "rpm_cmd",
        "Js",
        "Js_inv",
        "omg_fb",
        "ial_fb",
        "ibe_fb",
        "psi_mu_al_fb",
        "psi_mu_be_fb",
        "rotor_flux_cmd",
        "omg_ctrl_err",
        "speed_ctrl_err",
        "iMs",
        "iTs",
        "theta_M",
        "cosT",
        "sinT",
        "omega_syn",
        "omega_sl",
        "uMs_cmd",
        "uTs_cmd",
        "iMs_cmd",
        "iTs_cmd",
        "pi_speed",
        "pi_iMs",
        "pi_iTs",
    )

    def __init__(self):
        self.vc_count = 0

        self.timebase = 0.0

        self.ual = 0.0
        self.ube = 0.0

        self.rs = 0.0
        self.rreq = 0.0
        self.Lsigma = 0.0
        self.alpha = 0.0
        self.Lmu = 0.0
        self.Lmu_inv = 0.0

        self.Tload = 0.0
        self.rpm_cmd = 0.0

        self.Js = 0.0
        self.Js_inv = 0.0

        self.omg_fb = 0.0
        self.ial_fb = 0.0
        self.ibe_fb = 0.0
        self.psi_mu_al_fb = 0.0
        self.psi_mu_be_fb = 0.0

        self.rotor_flux_cmd = 0.0

        self.omg_ctrl_err = 0.0
        self.speed_ctrl_err = 0.0

        self.iMs = 0.0
        self.iTs = 0.0

        self.theta_M = 0.0
        self.cosT = 0.0
        self.sinT = 0.0

        self.omega_syn = 0.0
        self.omega_sl = 0.0

        self.uMs_cmd = 0.0
        self.uTs_cmd = 0.0
        self.iMs_cmd = 0.0
        self.iTs_cmd = 0.0

        self.pi_speed = PI_REG()
        self.pi_iMs = PI_REG()
        self.pi_iTs = PI_REG()

    def CTRL_INIT(self, IM):
        # int i=0,j=0
        self.cnt = 0
        self.vc_count = 0

        self.timebase = 0.0

        # Parameter (including speed) Adaptation */
        self.rs = IM.rs
        self.rreq = IM.rreq
        self.Lsigma = IM.Lsigma
        self.alpha = IM.alpha
        self.Lmu = IM.Lmu
        self.Lmu_inv = 1.0 / IM.Lmu
        self.Js = IM.Js
        self.Js_inv = 1.0 / IM.Js

        self.ual = 0.0
        self.ube = 0.0

        self.rpm_cmd = 0.0
        self.rotor_flux_cmd = 0.5

        self.omg_ctrl_err = 0.0
        self.speed_ctrl_err = 0.0

        self.omg_fb = 0.0
        self.ial_fb = 0.0
        self.ibe_fb = 0.0
        self.psi_mu_al_fb = 0.0
        self.psi_mu_be_fb = 0.0

        self.theta_M = 0.0
        self.cosT = 0.0
        self.sinT = 1

        self.uMs_cmd = 0.0
        self.uTs_cmd = 0.0
        self.iMs_cmd = 0.0
        self.iTs_cmd = 0.0

        self.omega_syn = 0.0
        self.omega_sl = 0.0

        # ver. IEMDC
        self.pi_speed.Kp = 10
        self.pi_speed.Ti = 0.8
        self.pi_speed.Ki = (
            (self.pi_speed.Kp * 4.77)
            / self.pi_speed.Ti
            * (TS * VC_LOOP_CEILING * DOWN_FREQ_EXE_INVERSE)
        )
        self.pi_speed.i_state = 0.0
        self.pi_speed.i_limit = 8

        print("Kp_omg={0}, Ki_omg={1}\n".format(self.pi_speed.Kp, self.pi_speed.Ki))

        self.pi_iMs.Kp = 15  # cutoff frequency of 1530 rad/s
        self.pi_iMs.Ti = 0.08
        self.pi_iMs.Ki = self.pi_iMs.Kp / self.pi_iMs.Ti * TS  # =0.025
        self.pi_iMs.i_state = 0.0
        self.pi_iMs.i_limit = 350  # 350.0  # unit: Volt

        self.pi_iTs.Kp = 15
        self.pi_iTs.Ti = 0.08
        self.pi_iTs.Ki = self.pi_iTs.Kp / self.pi_iTs.Ti * TS
        # self.pi_iTs.Ki = 0.0
        self.pi_iTs.i_state = 0.0
        self.pi_iTs.i_limit = 650  # unit: Volt, 350V->max 1300rpm

        print("Kp_cur={0}, Ki_cur={1}\n".format(self.pi_iMs.Kp, self.pi_iMs.Ki))

    def PI(self, r, err):
        # define I_STATE r->i_state
        # define I_LIMIT r->i_limit
        # double output;
        r.i_state += err * r.Ki  # 积分
        if r.i_state > r.i_limit:  # 添加积分饱和特性
            r.i_state = r.i_limit
        elif r.i_state < -r.i_limit:
            r.i_state = -r.i_limit

        output = r.i_state + err * r.Kp

        if output > r.i_limit:
            output = r.i_limit
        elif output < -r.i_limit:
            output = -r.i_limit
        return output

    def cmd_fast_speed_reversal(self, instant, interval, rpm_cmd, IM):
        if self.timebase > instant + 2 * interval:
            IM.rpm_cmd = rpm_cmd
        elif self.timebase > instant + interval:
            IM.rpm_cmd = -rpm_cmd
        elif self.timebase > instant:
            IM.rpm_cmd = rpm_cmd

    def control(self, speed_cmd, speed_cmd_dot, OB, IM):
        # OPEN LOOP CONTROL
        if CONTROL_STRATEGY == VVVF_CONTROL:
            self.VF_RATIO = 18  # 18.0 # 8 ~ 18 shows saturated phenomenon
            self.freq = 2  # 0.15 ~ 0.5 ~ 2 （0.1时电压李萨茹就变成一个圆了）
            self.volt = self.VF_RATIO * self.freq
            self.ual = self.volt * np.cos(2 * np.pi * self.freq * self.timebase)
            self.ube = self.volt * np.sin(2 * np.pi * self.freq * self.timebase)
            return

        # Input 1 is feedback: estimated speed or measured speed
        if SENSORLESS_CONTROL:
            self.omg_fb = OB.tajima.omg
            # #if OBSERVER_APPLIED == TAJIMA96
            #     CTRL.omega_syn = ob.tajima.omega_syn;
            #     CTRL.omega_sl  = ob.tajima.omega_sl;
            # #endif
        else:
            self.omg_fb = OB.im.omg

        # Input 2 is feedback: measured current
        self.ial_fb = OB.im.is_curr[0]
        self.ibe_fb = OB.im.is_curr[1]
        # Input 3 differs for DFOC and IFOC
        if CONTROL_STRATEGY == DFOC:
            # DFOC: estimated flux components in alpha-beta frame
            self.psi_mu_al_fb = OB.psi_mu_al[0]
            self.psi_mu_be_fb = OB.psi_mu_al[0]
        elif CONTROL_STRATEGY == IFOC:
            # IFOC: estimated rotor resistance
            self.rreq = OB.rreq

        # Flux (linkage) command
        self.rotor_flux_cmd = (
            1.5  # f(speed, dc bus voltage, last torque current command)
        )
        # CTRL.rotor_flux_cmd = 3;
        # 1. speed is compared with the base speed to decide flux weakening or not
        # 2. dc bus voltage is required for certain application
        # 3. last torque current command is required for loss minimization

        # M-axis current command
        self.iMs_cmd = (
            self.rotor_flux_cmd * self.Lmu_inv
            + M1 * OMG1 * np.cos(OMG1 * self.timebase) / self.rreq
        )
        # printf("%g, %g, %g\n", CTRL.Lmu_inv, CTRL.iMs_cmd, CTRL.iTs_cmd);

        # T-axis current command
        self.vc_count += 1
        if self.vc_count == VC_LOOP_CEILING * DOWN_FREQ_EXE_INVERSE:
            # if(True):
            self.vc_count = 0
            self.omg_ctrl_err = self.omg_fb - speed_cmd * IM.RPM_2_RAD_PER_SEC
            self.iTs_cmd = -self.PI(self.pi_speed, self.omg_ctrl_err)

            self.speed_ctrl_err = self.omg_ctrl_err * IM.RAD_PER_SEC_2_RPM

        if CONTROL_STRATEGY == DFOC:
            # feedback field orientation
            modulus = np.sqrt(
                self.psi_mu_al_fb * self.psi_mu_al_fb
                + self.psi_mu_be_fb * self.psi_mu_be_fb
            )
            if modulus < 1e-3:
                self.cosT = 1
                self.sinT = 0
            else:
                self.cosT = self.psi_mu_al_fb  # modulus
                self.sinT = self.psi_mu_be_fb  # modulus
        elif CONTROL_STRATEGY == IFOC:
            # Feed-forward field orientation
            self.theta_M += TS * self.omega_syn

            if self.theta_M > np.pi:
                self.theta_M -= 2 * np.pi
            elif self.theta_M < -np.pi:
                self.theta_M += 2 * np.pi  # 反转！

            if VOLTAGE_CURRENT_DECOUPLING_CIRCUIT:
                self.omega_sl = self.rreq * self.iTs / self.rotor_flux_cmd
            else:
                self.omega_sl = self.rreq * self.iTs_cmd / self.rotor_flux_cmd

            self.omega_syn = self.omg_fb + self.omega_sl

            self.cosT = np.cos(self.theta_M)
            self.sinT = np.sin(self.theta_M)

        # Measured current in M-T frame
        self.iMs = mc.AB2M(self.ial_fb, self.ibe_fb, self.cosT, self.sinT)
        self.iTs = mc.AB2T(self.ial_fb, self.ibe_fb, self.cosT, self.sinT)

        # Voltage command in M-T frame
        vM = -self.PI(self.pi_iMs, self.iMs - self.iMs_cmd)
        vT = -self.PI(self.pi_iTs, self.iTs - self.iTs_cmd)

        # Current loop decoupling (skipped, see Chen.Huang-Stable)
        #   # Steady state dynamics based decoupling circuits for current regulation
        if VOLTAGE_CURRENT_DECOUPLING_CIRCUIT == True:
            self.uMs_cmd = vM + (self.Lsigma) * (
                -self.omega_syn * self.iTs
            )  # Telford03/04
            self.uTs_cmd = vT + self.omega_syn * (
                self.rotor_flux_cmd + self.Lsigma * self.iMs
            )  # 这个行，但是无速度运行时，会导致M轴电流在转速暂态高频震荡。
            # CTRL.uMs_cmd = vM;
            # CTRL.uTs_cmd = vT;
        else:
            self.uMs_cmd = vM
            self.uTs_cmd = vT

        # Voltage command in alpha-beta frame
        self.ual = mc.MT2A(self.uMs_cmd, self.uTs_cmd, self.cosT, self.sinT)
        self.ube = mc.MT2B(self.uMs_cmd, self.uTs_cmd, self.cosT, self.sinT)
