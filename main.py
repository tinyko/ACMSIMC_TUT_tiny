from ACM import ACIM
from CTRL import CTRL0,s_curve
from Observer import Observer
from Utils import *
from FileIO import FileIO
import Macros as mc
import ACMPlot

#Debug mode raise all runtime warning
np.seterr(all='raise')

def main():
    IM=ACIM()
    OB=Observer()
    CTRL=CTRL0()
    IM.Machine_init()
    OB.acm_init(IM)
    OB.ob_init()
    CTRL.CTRL_INIT(IM)
    sc=s_curve()

    IM_items=vars(IM)
    OB_items=vars(OB)
    CTRL_items=vars(CTRL)

    # print("{0}\n{1}\n{2}\n".format(IM_items,OB_items,CTRL_items))

    dfe=0 # dfe for down frequency execution

    f=open(r'algorithm.dat','w')
    fio=FileIO()
    fio.write_header_to_file(f)

    for _ in range(NUMBER_OF_LINES):

        # Command and Load Torque */
        # cmd_fast_speed_reversal(CTRL.timebase, 5, 5, 1500); // timebase, instant, interval, rpm_cmd

        # CTRL.cmd_fast_speed_reversal(5, 5, 100, IM) # timebase, instant, interval, rpm_cmd
        sc.speed_ref(CTRL.timebase,IM)
        # ACM.Tload = 5 * sign(ACM.rpm); 
        # if CTRL.timebase>= 5.00275 :
        #     print('no')
        
        IM.Tload = 5 * np.sign(IM.rpm) # No-load test
        #print(IM.rpm)
        # ACM.Tload = ACM.Tem; // Blocked-rotor test

        # Simulated ACM */
        if(IM.machine_simulation(CTRL)): 
            print("Break the loop.\n")
            break
        dfe+=1
        if(dfe==DOWN_FREQ_EXE):
            dfe = 0

            # Time */
            CTRL.timebase += TS

            OB.measurement(IM,CTRL)

            OB.observation(CTRL)

            fio.write_data_to_file(f,IM,CTRL)
             
            CTRL.control(IM.rpm_cmd, 0 ,OB,IM)
        
        IM.inverter_model(CTRL)
    f.close()

    ACMPlot.draw_trend()

if __name__ == "__main__" :
    main()