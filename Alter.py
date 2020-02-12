import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
np.seterr(all='raise')
INDUCTION_MACHINE=1
SYNCHRONOUS_MACHINE=2
MACHINE_TYPE=INDUCTION_MACHINE

VVVF_CONTROL=0
IFOC=1
DFOC=2
CONTROL_STRATEGY=IFOC

TAJIMA96=0
TOTALLY_ADAPTIVE_OBSERVER=1 # Flux-MRAS
OBSERVER_APPLIED=1

SENSORLESS_CONTROL=False
VOLTAGE_CURRENT_DECOUPLING_CIRCUIT=True
SATURATED_MAGNETIC_CIRCUIT=True
INVERTER_NONLINEARITY=True

RUN_TIME=15.0 #in seconds
MACHINE_TS =1.25e-4
MACHINE_TS_INVERSE =8000
NUMBER_OF_LINES =int((RUN_TIME/MACHINE_TS))
DOWN_FREQ_EXE =2
DOWN_FREQ_EXE_INVERSE =0.5
TS=(MACHINE_TS*DOWN_FREQ_EXE) #2.5e-4 
TS_INVERSE =(MACHINE_TS_INVERSE*DOWN_FREQ_EXE_INVERSE) # 4000
DOWN_SAMPLE =1 # 5 # 10

selfSIMC_DEBUG=False

IZ_STEP=0.01
IZ_STEP_INV=100

one_over_six = 1.0/6.0 

# Macro for Park transformation*/
AB2M=lambda A, B, COS, SIN : ( (A)*COS  + (B)*SIN )
AB2T=lambda A, B, COS, SIN : ( (A)*-SIN + (B)*COS ) 
MT2A=lambda M, T, COS, SIN : ( (M)*COS - (T)*SIN )
MT2B=lambda M, T, COS, SIN : ( (M)*SIN + (T)*COS )

isNumber = lambda x : isinstance(x,(np.int,np.int64,np.float,np.float64))

# Macro for Power-invariant inverse Clarke transformation
AB2U=lambda A, B : ( 0.816496580927726 * ( A ) )
AB2V=lambda A, B : ( 0.816496580927726 * ( A*-0.5 + B*0.8660254037844387 ) )
AB2W=lambda A, B : ( 0.816496580927726 * ( A*-0.5 + B*-0.8660254037844385 ) )

# General Constants */
lTpi=0.15915494309189535 # 1/(2*pi)
TWO_PI_OVER_3=2.0943951023931953
SIN_2PI_SLASH_3=0.86602540378443871 # sin(2*pi/3)
SIN_DASH_2PI_SLASH_3=-0.86602540378443871 # sin(-2*pi/3)
SQRT_2_SLASH_3=0.81649658092772603 # sqrt(2.0/3.0)
#  PI_D 3.1415926535897932384626433832795 #  */
M_PI_OVER_180=0.017453292519943295

# 补充的宏，为了实现实验/仿真代码大一统
PHASE_NUMBER=3

_Vce0 =1.8 # V
_Vd0  =1.3 # V
_Udc  =300 # V
_Toff =0.32e-6 # sec
_Ton  =0.15e-6 # sec
_Tdead=3.30e-6 # sec
_Tcomp=0.0*(_Tdead+_Ton-_Toff) #3.13e-6  #8e-6  # 过补偿 #3.13e-6  # 只补偿死
Rce=0.04958
Rdiode=0.05618

VC_LOOP_CEILING=4

M1=0
OMG1=2

class ACIM():
    def __init__(self):
        self.ial=0.0 #Current Alpha Component
        self.ibe=0.0 #Current Beta Component
        self.psi_al=0.0 #
        self.psi_be=0.0 #
        self.ual=0.0 #Voltage Alpha Component
        self.ube=0.0 #Voltage Beta Component

        self.x=np.zeros(10) #Motor Status
        self.rpm=0.0 #revolutions per minute
        self.rpm_cmd=0.0 #speed command
        self.rpm_deriv_cmd=0.0 #derivative of speed command 
        self.Tload=0.0 #Load Torque
        self.Tem=0.0 #Electrical Torque

        self.Js=0.0 #moment of inertia
        self.npp=0.0 #number of pole pairs
        self.mu_m=0.0 #=npp/Js
        self.Ts=0.0 #sampling time

        self.Lsigma=0.0 #Leakage induction
        self.rs=0.0 #stator resistance
        self.rreq=0.0 #equivalent rotor resistance
        self.Lmu=0.0 #magnetic induction
        self.Lmu_inv=0.0 #inverse of Lmu
        self.alpha=0.0 #inverse of rotor time constant = Lmu/Lr

        self.Lm=0.0 #mutual induction
        self.rr=0.0 #rotor resistance
        self.Lls=0.0 #stator leakage induction
        self.Llr=0.0 #rotor leakage induction
        self.Lm_slash_Lr=0.0
        self.Lr_slash_Lm=0.0

        self.LSigmal=0.0

        self.izq=0.0
        self.izd=0.0
        self.iz=0.0
        
        self.psimq=0.0
        self.psimd=0.0
        self.psim=0.0
        self.im=0.0

        self.iqs=0.0
        self.ids=0.0
        self.iqr=0.0
        self.idr=0.0

        self.ual_c_dist=0.0
        self.ube_c_dist=0.0
        self.dist_al=0.0
        self.dist_be=0.0
        
        self.RAD_PER_SEC_2_RPM=0.0
        self.RPM_2_RAD_PER_SEC=0.0
        
    def Machine_init(self):
        self.db=[]
        for i in range(5):
            self.x[i] = 0.0
        self.rpm = 0.0 
        self.rpm_cmd = 0.0 
        self.rpm_deriv_cmd = 0.0 
        self.Tload = 0.0 
        self.Tem = 0.0 

        self.Lmu    = 0.4482 
        # Those parameters are introduced since branch saturation
        self.Lls = 0.0126 
        self.Llr = 0.0126 
        self.Lm  = 0.5*(self.Lmu + np.sqrt(self.Lmu*self.Lmu + 4*self.Llr*self.Lmu)) 
        self.Lm_slash_Lr = self.Lm/(self.Lm+self.Llr) 
        self.Lr_slash_Lm = (self.Lm+self.Llr)/self.Lm   # i_dreq = self.idr*self.Lr_slash_Lm
        self.LSigmal = 1.0 / (1.0 / self.Lls + 1.0 / self.Llr) 
        self.Lsigma = self.Lm + self.Lls - self.Lmu   # = Ls * (1.0 - Lm*Lm/Ls/Lr) 
        print("Validate: {0} = {1} ?\n".format(self.Lsigma, (self.Lm+self.Lls) * (1.0 - self.Lm*self.Lm/(self.Lm+self.Lls)/(self.Lm+self.Llr))) ) 

        self.rreq   = 1.69 
        self.rs     = 3.04 
        self.rr = self.rreq * self.Lr_slash_Lm*self.Lr_slash_Lm 

        self.alpha  = self.rreq / (self.Lmu) 
        self.Lmu_inv= 1.0/self.Lmu 

        self.Js = 0.0636   # Awaya92 using im.omg
        self.npp = 2 
        self.mu_m = self.npp/self.Js 

        self.Ts  = MACHINE_TS 

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
        
        self.satLUT=np.array(pd.read_csv(r"33.txt",header=None))
        self.RAD_PER_SEC_2_RPM=(60.0/(2*np.pi*self.npp))
        self.RPM_2_RAD_PER_SEC=((2*np.pi*self.npp)/60.0)
        
    def sat_lookup(self,iz,satLUT):
        if iz >= (444000)*IZ_STEP :
            idx = 444000
            return (iz - idx*IZ_STEP) * IZ_STEP_INV * (satLUT[idx+1][0]-satLUT[idx][0]) + satLUT[idx][0] # 一阶插值？
        # return satLUT[(int)(iz/IZ_STEP+0.5)]  # This yields bad results with obvious oscillation (voltages and currents)
        idx = int(iz*IZ_STEP_INV)
        return (iz - idx*IZ_STEP) * IZ_STEP_INV * (satLUT[idx+1][0]-satLUT[idx][0]) + satLUT[idx][0] # 一阶插值？
        
    def collectCurrents(self):
        # Generalised Current by Therrien2013
        self.izq = self.x[1]/self.Lls + self.x[3]/self.Llr 
        self.izd = self.x[0]/self.Lls + self.x[2]/self.Llr
        
        self.iz = np.sqrt(self.izd*self.izd + self.izq*self.izq)

        if self.iz>1e-8:
            if SATURATED_MAGNETIC_CIRCUIT:
                self.psim = self.sat_lookup(self.iz, self.satLUT)
                self.im = self.iz - self.psim/self.LSigmal 
                self.Lm = self.psim/self.im 
                self.Lmu = self.Lm*self.Lm/(self.Lm+self.Llr) 
                self.Lmu_inv = 1.0/self.Lmu 
                self.alpha = self.rr/(self.Lm+self.Llr) 
                self.rreq = self.Lmu*self.alpha 
                self.Lsigma = (self.Lls+self.Lm) - self.Lmu 
                self.Lm_slash_Lr = self.Lm/(self.Lm+self.Llr) 
                self.Lr_slash_Lm = (self.Lm+self.Llr)/self.Lm 
            else :
                self.psim = 1.0/(1.0/self.Lm+1.0/self.Lls+1.0/self.Llr)*self.iz 
            self.psimq = self.psim/self.iz*self.izq 
            self.psimd = self.psim/self.iz*self.izd 
        else :
            if selfSIMC_DEBUG :
                print("how to handle zero iz?\n")
            self.psimq = 0.0 
            self.psimd = 0.0 
        self.iqs = (self.x[1] - self.psimq) / self.Lls 
        self.ids = (self.x[0] - self.psimd) / self.Lls 
        self.iqr = (self.x[3] - self.psimq) / self.Llr 
        self.idr = (self.x[2] - self.psimd) / self.Llr     
        # # Direct compute is ir from psis psir */
        # self.iqs = (self.x[1] - (self.Lm/(self.Lm+self.Llr))*self.x[3])/self.Lsigma 
        # self.ids = (self.x[0] - (self.Lm/(self.Lm+self.Llr))*self.x[2])/self.Lsigma 
        # self.iqr = (self.x[3] - (self.Lm/(self.Lm+self.Lls))*self.x[1])/(self.Lm+self.Llr-self.Lm*self.Lm/(self.Lm+self.Lls)) 
        # self.idr = (self.x[2] - (self.Lm/(self.Lm+self.Lls))*self.x[0])/(self.Lm+self.Llr-self.Lm*self.Lm/(self.Lm+self.Lls)) 
    
    def rK5_satDynamics(self,t, x, fx):
        # argument t is omitted*/

        # STEP ZERO: collect all the currents: is, ir, iz */
        self.collectCurrents()

        # STEP ONE: Inverter Nonlinearity Considered */
        # not CTRL.ual
        # not CTRL.ube

        # STEP TWO: electromagnetic model with full flux states in alpha-beta frame */
        fx[1] = self.ube - self.rs*self.iqs  # 这里是反过来的！Q轴在前面！
        fx[0] = self.ual - self.rs*self.ids 
        fx[3] =          - self.rr*self.iqr + x[4]*x[2]  # 这里是反过来的！Q轴在前面！
        fx[2] =          - self.rr*self.idr - x[4]*x[3] 

        # STEP THREE: mechanical model */
        # self.Tem = self.npp*(self.Lm/(self.Lm+self.Llr))*(self.iqs*self.x[2]-self.ids*self.x[3])  # this is not better 
        self.Tem = self.npp*(self.iqs*x[0]-self.ids*x[1]) 
        fx[4] = (self.Tem - self.Tload)*self.mu_m 
    
    def rK555_Sat(self,t, x, hs):
        k1=np.zeros(5)
        k2=np.zeros(5)
        k3=np.zeros(5)
        k4=np.zeros(5)
        xk=np.zeros(5)
        fx=np.zeros(5)

        self.rK5_satDynamics(t, x, fx)  # timer.t,
        for i in range(5):    
            k1[i] = fx[i] * hs 
            xk[i] = x[i] + k1[i]*0.5 

        self.rK5_satDynamics(t, xk, fx)  # timer.t+hs/2., 
        for i in range(5):        
            k2[i] = fx[i] * hs 
            xk[i] = x[i] + k2[i]*0.5 

        self.rK5_satDynamics(t, xk, fx)  # timer.t+hs/2., 
        for i in range(5):      
            k3[i] = fx[i] * hs 
            xk[i] = x[i] + k3[i] 

        self.rK5_satDynamics(t, xk, fx)  # timer.t+hs, 
        for i in range(5):      
            k4[i] = fx[i] * hs 
            x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])*one_over_six 
        self.collectCurrents()

    def rK5_dynamics(self,t, x, fx):
        if MACHINE_TYPE == INDUCTION_MACHINE:
            # electromagnetic model
            fx[2] = self.rreq*x[0] - self.alpha*x[2] - x[4]*x[3]  # flux-alpha
            fx[3] = self.rreq*x[1] - self.alpha*x[3] + x[4]*x[2]  # flux-beta
            fx[0] = (self.ual - self.rs*x[0] - fx[2])/self.Lsigma  # current-alpha
            fx[1] = (self.ube - self.rs*x[1] - fx[3])/self.Lsigma  # current-beta

            # mechanical model
            self.Tem = self.npp*(x[1]*x[2]-x[0]*x[3]) 
            fx[4] = (self.Tem - self.Tload)*self.mu_m  # elec. angular rotor speed
            #fx[5] = x[4]                            # elec. angular rotor position

    def rK555_Lin(self,t, x, hs):
        if MACHINE_TYPE == INDUCTION_MACHINE:
            NUMBER_OF_STATES=6
        NS=NUMBER_OF_STATES
        k1=np.zeros(NS)
        k2=np.zeros(NS)
        k3=np.zeros(NS)
        k4=np.zeros(NS)
        xk=np.zeros(NS)
        fx=np.zeros(NS)
        
        self.rK5_dynamics(t, x, fx)  # timer.t,
        for i in range(NS):        
            k1[i] = fx[i] * hs 
            xk[i] = x[i] + k1[i]*0.5 

        self.rK5_dynamics(t, xk, fx)  # timer.t+hs/2., 
        for i in range(NS):        
            k2[i] = fx[i] * hs 
            xk[i] = x[i] + k2[i]*0.5 

        self.rK5_dynamics(t, xk, fx)  # timer.t+hs/2., 
        for i in range(NS):        
            k3[i] = fx[i] * hs 
            xk[i] = x[i] + k3[i] 

        self.rK5_dynamics(t, xk, fx)  # timer.t+hs, 
        for i in range(NS):        
            k4[i] = fx[i] * hs 
            x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0 
            
    def machine_simulation(self):

        # API for explicit access
        if MACHINE_TYPE == INDUCTION_MACHINE:
            # self.rK555_Lin(CTRL.timebase, self.x, self.Ts) 
            # self.ial    = self.x[0]  # rK555_Lin
            # self.ibe    = self.x[1]  # rK555_Lin
            # self.psi_al = self.x[2]  # rK555_Lin
            # self.psi_be = self.x[3]  # rK555_Lin

           
            self.rK555_Sat(CTRL.timebase, self.x, self.Ts) 
            self.ial    = self.ids  # rK555_Sat
            self.ibe    = self.iqs  # rK555_Sat
            self.psi_al = self.x[2]*self.Lm_slash_Lr  # rK555_Sat
            self.psi_be = self.x[3]*self.Lm_slash_Lr  # rK555_Sat

            self.rpm    = self.x[4] * 60 / (2 * np.pi * self.npp) 


        if(self.rpm!=np.nan):
            return False
        else :
            print("self.rpm is {0}\n", self.rpm)
            return True       

    def InverterNonlinearity_SKSul96(self, ual, ube, ial, ibe):
        TM = _Toff - _Ton - _Tdead + _Tcomp # Sul1996
        Udist = (_Udc*TM*TS_INVERSE - _Vce0 - _Vd0) / 6.0 # Udist = (_Udc*TM/1e-4 - _Vce0 - _Vd0) / 6.0 
        # Udist = (_Udc*TM*TS_INVERSE) / 6.0 
        # Udist = 0.0 

        ia = SQRT_2_SLASH_3 * (       ial                              )
        ib = SQRT_2_SLASH_3 * (-0.5 * ial - SIN_DASH_2PI_SLASH_3 * ibe )
        ic = SQRT_2_SLASH_3 * (-0.5 * ial - SIN_2PI_SLASH_3      * ibe )

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
        self.ual_c_dist = ual + np.sqrt(1.5)*Udist*(2*np.sign(ia) - np.sign(ib) - np.sign(ic)) - 0.5*(Rce+Rdiode)*ial
        self.ube_c_dist = ube + 3/np.sqrt(2)*Udist*(                np.sign(ib) - np.sign(ic)) - 0.5*(Rce+Rdiode)*ibe 

    def inverter_model(self):

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

        # 冗余变量赋值

def measurement():
    OB.im.us_curr[0] = CTRL.ual
    OB.im.us_curr[1] = CTRL.ube
    OB.im.us_prev[0] = OB.im.us_curr[0]
    OB.im.us_prev[1] = OB.im.us_curr[1]

    if MACHINE_TYPE == INDUCTION_MACHINE:
        OB.im.is_curr[0] = IM.ial
        OB.im.is_curr[1] = IM.ibe
        OB.im.omg = IM.x[4]
        OB.im.theta_r = IM.x[5]

class PI_REG():
    def __init__(self):
        self.Kp=0.0
        self.Ti=0.0
        self.Ki=0.0 # Ki = Kp/Ti*TS
        self.i_state=0.0
        self.i_limit=0.0

class Tajima():
    def __init__(self):
        self.K_PEM=0.0
        self.omega_syn=0.0
        self.e_M=0.0
        self.omega_sl=0.0
        self.omg=0.0

class InductionMachine():
    def __init__(self):
        self.zeros=np.zeros(2)
        self.u_s=self.zeros
        self.i_s=self.zeros
        self.us_curr=self.zeros
        self.is_curr=self.zeros
        self.us_prev=self.zeros
        self.is_prev=self.zeros

        self.Js=0.0
        self.Js_inv=0.0

        self.rs=0.0
        self.rreq=0.0

        self.alpha=0.0
        self.Lsigma=0.0
        self.Lsigma_inv=0.0
        self.Lmu=0.0
        self.Lmu_inv=0.0

        self.npp=0.0
        self.omg=0.0
        self.theta_r=0.0
        self.theta_d=0.0

class Observer():
    def __init__(self):
        self.Js=0.0
        self.Js_inv=0.0

        self.rs=0.0
        self.rreq=0.0

        self.alpha=0.0
        self.Lsigma=0.0
        self.Lsigma_inv=0.0
        self.Lmu=0.0
        self.Lmu_inv=0.0

        self.npp=0.0
        self.omg=0.0

        self.psi_mu_al=0.0
        self.psi_mu_be=0.0

        self.tajima=Tajima()
        self.im=InductionMachine()
        
    def acm_init(self):
        self.im.u_s = self.im.zeros
        self.im.i_s = self.im.zeros
        self.im.us_curr = self.im.zeros
        self.im.is_curr = self.im.zeros
        self.im.us_prev = self.im.zeros
        self.im.is_prev = self.im.zeros      

        self.im.Js = IM.Js
        self.im.Js_inv = 1./self.im.Js
        self.im.rs = IM.rs
        self.im.rreq = IM.rreq
        self.im.alpha = IM.alpha
        self.im.Lsigma = IM.Lsigma
        self.im.Lsigma_inv = 1/self.im.Lsigma
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
    def observation(self):
        rotor_flux_cmd = CTRL.rotor_flux_cmd
        iMs = CTRL.iMs
        iTs = CTRL.iTs
        uMs_cmd  = CTRL.uMs_cmd
        uTs_cmd  = CTRL.uTs_cmd

        # Speed estimation: Tajima1996
        self.tajima.K_PEM = 0.1 # 0.1 ~ 2
        self.tajima.omega_syn  = (uTs_cmd - self.rs * iTs) / (rotor_flux_cmd + self.Lsigma * iMs)
        self.tajima.e_M        = uMs_cmd - self.rs * iMs + self.tajima.omega_syn * self.Lsigma * iTs
        self.tajima.omega_syn -= self.tajima.K_PEM * self.tajima.e_M

        self.tajima.omega_sl   = self.rreq * iTs / rotor_flux_cmd
        self.tajima.omg = self.tajima.omega_syn - self.tajima.omega_sl # Instantaneous Velocity Computation

        # Flux estimation 1: Voltage model (implemented by shitty Euler method for now)
        
        deriv_psi_s_al = self.im.us_curr[0] - self.rs*self.im.is_curr[0]
        deriv_psi_s_be = self.im.us_curr[1] - self.rs*self.im.is_curr[1]
        psi_s_al=psi_s_be=0
        
        self.psi_s_al += TS*deriv_psi_s_al
        self.psi_s_be += TS*deriv_psi_s_be
        self.psi_mu_al = self.psi_s_al - self.Lsigma*self.im.i_s[0]
        self.psi_mu_be = self.psi_s_be - self.Lsigma*self.im.i_s[1]

        # Flux estimation 2: Current model (this is a bad flux estimator)  

class CTRL0():
    def __init__(self):
        self.vc_count = 0
        
        self.timebase=0.0 

        self.ual=0.0  
        self.ube=0.0  

        self.rs=0.0  
        self.rreq=0.0  
        self.Lsigma=0.0  
        self.alpha=0.0  
        self.Lmu=0.0  
        self.Lmu_inv=0.0  

        self.Tload=0.0  
        self.rpm_cmd=0.0  

        self.Js=0.0  
        self.Js_inv=0.0  

        self.omg_fb=0.0  
        self.ial_fb=0.0  
        self.ibe_fb=0.0  
        self.psi_mu_al_fb=0.0  
        self.psi_mu_be_fb=0.0  

        self.rotor_flux_cmd=0.0  

        self.omg_ctrl_err=0.0  
        self.speed_ctrl_err=0.0  

        self.iMs=0.0  
        self.iTs=0.0  

        self.theta_M=0.0  
        self.cosT=0.0  
        self.sinT=0.0  

        self.omega_syn=0.0  
        self.omega_sl=0.0  

        self.uMs_cmd=0.0  
        self.uTs_cmd=0.0  
        self.iMs_cmd=0.0  
        self.iTs_cmd=0.0  

        self.pi_speed=PI_REG()
        self.pi_iMs=PI_REG()
        self.pi_iTs=PI_REG()
    def CTRL_INIT(self):
        #int i=0,j=0 
        self.cnt=0
        self.vc_count = 0
        
        self.timebase = 0.0 

        # Parameter (including speed) Adaptation */ 
        self.rs     = IM.rs 
        self.rreq   = IM.rreq 
        self.Lsigma = IM.Lsigma 
        self.alpha  = IM.alpha 
        self.Lmu    = IM.Lmu 
        self.Lmu_inv = 1.0/IM.Lmu 
        self.Js     = IM.Js 
        self.Js_inv = 1.0/IM.Js     

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
        self.pi_speed.Kp = 0.5  
        self.pi_speed.Ti = 5 
        self.pi_speed.Ki = (self.pi_speed.Kp*4.77) / self.pi_speed.Ti * (TS*VC_LOOP_CEILING*DOWN_FREQ_EXE_INVERSE) 
        self.pi_speed.i_state = 0.0 
        self.pi_speed.i_limit = 8 

        print("Kp_omg={0}, Ki_omg={1}\n".format(self.pi_speed.Kp, self.pi_speed.Ki)) 

        self.pi_iMs.Kp = 15  # cutoff frequency of 1530 rad/s
        self.pi_iMs.Ti = 0.08 
        self.pi_iMs.Ki = self.pi_iMs.Kp/self.pi_iMs.Ti*TS  # =0.025
        self.pi_iMs.i_state = 0.0 
        self.pi_iMs.i_limit = 350  #350.0  # unit: Volt

        self.pi_iTs.Kp = 15 
        self.pi_iTs.Ti = 0.08 
        self.pi_iTs.Ki = self.pi_iTs.Kp/self.pi_iTs.Ti*TS 
        self.pi_iTs.i_state = 0.0 
        self.pi_iTs.i_limit = 650  # unit: Volt, 350V->max 1300rpm

        print("Kp_cur={0}, Ki_cur={1}\n".format(self.pi_iMs.Kp, self.pi_iMs.Ki))
        
    def PI(self, r, err):
        #define I_STATE r->i_state
        #define I_LIMIT r->i_limit
        #double output;
        r.i_state += err * r.Ki    # 积分
        if( r.i_state > r.i_limit):    # 添加积分饱和特性
            r.i_state = r.i_limit 
        elif( r.i_state < -r.i_limit):
            r.i_state = -r.i_limit

        output = r.i_state + err * r.Kp

        if(output > r.i_limit):
            output = r.i_limit
        elif(output < -r.i_limit):
            output = -r.i_limit
        return output
    
    def control(self, speed_cmd, speed_cmd_dot):
        # OPEN LOOP CONTROL
        if CONTROL_STRATEGY == VVVF_CONTROL :
            self.VF_RATIO = 18 #18.0 # 8 ~ 18 shows saturated phenomenon
            self.freq = 2 # 0.15 ~ 0.5 ~ 2 （0.1时电压李萨茹就变成一个圆了）
            self.volt = self.VF_RATIO*self.freq
            self.ual = self.volt * np.cos(2* np.pi * self.freq * self.timebase)
            self.ube = self.volt * np.sin(2* np.pi * self.freq * self.timebase)
            return

        # Input 1 is feedback: estimated speed or measured speed
        if SENSORLESS_CONTROL :
            self.omg_fb = OB.tajima.omg
            # #if OBSERVER_APPLIED == TAJIMA96
            #     CTRL.omega_syn = ob.tajima.omega_syn;
            #     CTRL.omega_sl  = ob.tajima.omega_sl;
            # #endif
        else :
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
        self.rotor_flux_cmd = 1.5 # f(speed, dc bus voltage, last torque current command)
        # CTRL.rotor_flux_cmd = 3;
            # 1. speed is compared with the base speed to decide flux weakening or not
            # 2. dc bus voltage is required for certain application
            # 3. last torque current command is required for loss minimization

        # M-axis current command
        self.iMs_cmd = self.rotor_flux_cmd*self.Lmu_inv + M1 * OMG1 * np.cos(OMG1 * self.timebase) / self.rreq
        # printf("%g, %g, %g\n", CTRL.Lmu_inv, CTRL.iMs_cmd, CTRL.iTs_cmd);

        # T-axis current command
        self.vc_count+=1
        if(self.vc_count == VC_LOOP_CEILING * DOWN_FREQ_EXE_INVERSE):
        #if(True):
            self.vc_count = 0
            self.omg_ctrl_err = self.omg_fb - speed_cmd * IM.RPM_2_RAD_PER_SEC
            self.iTs_cmd = - self.PI(self.pi_speed, self.omg_ctrl_err)

            self.speed_ctrl_err = self.omg_ctrl_err * IM.RAD_PER_SEC_2_RPM

        if CONTROL_STRATEGY == DFOC:
            # feedback field orientation
            modulus = np.sqrt(self.psi_mu_al_fb*self.psi_mu_al_fb + self.psi_mu_be_fb*self.psi_mu_be_fb)
            if(modulus<1e-3):
                self.cosT = 1
                self.sinT = 0
            else:
                self.cosT = self.psi_mu_al_fb # modulus
                self.sinT = self.psi_mu_be_fb # modulus
        elif CONTROL_STRATEGY == IFOC:
            # Feed-forward field orientation
            self.theta_M += TS * self.omega_syn

            if(self.theta_M > np.pi):
                self.theta_M -= 2*np.pi
            elif(self.theta_M < -np.pi):
                self.theta_M += 2*np.pi # 反转！

            # if VOLTAGE_CURRENT_DECOUPLING_CIRCUIT:
            #     self.omega_sl = self.rreq*self.iTs / self.rotor_flux_cmd
            # else:
            self.omega_sl = self.rreq*self.iTs_cmd / self.rotor_flux_cmd

            self.omega_syn = self.omg_fb + self.omega_sl

            self.cosT = np.cos(self.theta_M)
            self.sinT = np.sin(self.theta_M)

        # Measured current in M-T frame
        self.iMs = AB2M(self.ial_fb, self.ibe_fb, self.cosT, self.sinT)
        self.iTs = AB2T(self.ial_fb, self.ibe_fb, self.cosT, self.sinT)

        # Voltage command in M-T frame
        vM = - self.PI(self.pi_iMs, self.iMs-self.iMs_cmd)
        vT = - self.PI(self.pi_iTs, self.iTs-self.iTs_cmd)

        # Current loop decoupling (skipped, see Chen.Huang-Stable)
        #   # Steady state dynamics based decoupling circuits for current regulation
        if VOLTAGE_CURRENT_DECOUPLING_CIRCUIT == True:
            self.uMs_cmd = vM + (self.Lsigma) * (-self.omega_syn * self.iTs) # Telford03/04
            self.uTs_cmd = vT + self.omega_syn * (self.rotor_flux_cmd + self.Lsigma*self.iMs) # 这个行，但是无速度运行时，会导致M轴电流在转速暂态高频震荡。
            # CTRL.uMs_cmd = vM;
            # CTRL.uTs_cmd = vT;
        else:
            self.uMs_cmd = vM
            self.uTs_cmd = vT

        # Voltage command in alpha-beta frame
        self.ual = MT2A(self.uMs_cmd, self.uTs_cmd, self.cosT, self.sinT)
        self.ube = MT2B(self.uMs_cmd, self.uTs_cmd, self.cosT, self.sinT)

def cmd_fast_speed_reversal(timebase, instant, interval, rpm_cmd):
    if(timebase > instant+2*interval):
        IM.rpm_cmd = rpm_cmd
    elif(timebase > instant+interval):
        IM.rpm_cmd = -rpm_cmd
    elif(timebase > instant):
        IM.rpm_cmd = rpm_cmd

class FileIO():
    def __init__(self):
        self.j=0
    def write_header_to_file(self,f):
        f.write("ACM.x[0],ACM.x[1],ACM.x[2],ACM.x[3],ACM.x[4],ACM.Tem,ACM.ual,ACM.ube,ACM.ial,ACM.ibe,ACM.rpm_cmd,e_omega\n")

        f2 = open(r"info.dat", "w")
        f2.write("TS,DOWN_SAMPLE\n")
        f2.write("{0}, {1}\n".format(TS, DOWN_SAMPLE))
        f2.close()
    def write_data_to_file(self,f):
        self.j+=1
        if(self.j == DOWN_SAMPLE):
            self.j=0
            f.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n".
                    format(IM.x[0],IM.x[1],IM.x[2],IM.x[3],IM.x[4],IM.Tem,IM.ual,IM.ube,IM.ial,IM.ibe,IM.rpm_cmd,IM.rpm_cmd-IM.x[4]*IM.RAD_PER_SEC_2_RPM))

IM=ACIM()
OB=Observer()
CTRL=CTRL0()
IM.Machine_init()
OB.acm_init()
OB.ob_init()
CTRL.CTRL_INIT()

dfe=0 # dfe for down frequency execution

f=open(r'algorithm.dat','w')
fio=FileIO()
fio.write_header_to_file(f)

for _ in range(NUMBER_OF_LINES):

    # Command and Load Torque */
    # cmd_fast_speed_reversal(CTRL.timebase, 5, 5, 1500); // timebase, instant, interval, rpm_cmd
    cmd_fast_speed_reversal(CTRL.timebase, 5, 5, 100) # timebase, instant, interval, rpm_cmd
    # ACM.Tload = 5 * sign(ACM.rpm); 
    # if CTRL.timebase>= 5.00275 :
    #     print('no')
    
    IM.Tload = 0 * np.sign(IM.rpm) # No-load test
    #print(IM.rpm)
    # ACM.Tload = ACM.Tem; // Blocked-rotor test

    # Simulated ACM */
    if(IM.machine_simulation()): 
        print("Break the loop.\n")
        break
    dfe+=1
    if(dfe==DOWN_FREQ_EXE):
        dfe = 0

        # Time */
        CTRL.timebase += TS

        measurement()

        OB.observation()

        fio.write_data_to_file(f)
        
        CTRL.control(IM.rpm_cmd, 0)
    

    IM.inverter_model()
f.close()