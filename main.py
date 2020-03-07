from ACM import ACIM
from CTRL import CTRL0, s_curve
from Observer import Observer
from Utils import np, NUMBER_OF_LINES, DOWN_FREQ_EXE, TS, RUN_TIME, MACHINE_TS_INVERSE
from FileIO import FileIO

# import Macros as mc
import ACMPlot2
import time

# import cProfile

# Debug mode raise all runtime warning
np.seterr(all="raise")


def main():
    start_time = time.time()
    IM = ACIM()
    IM2 = ACIM()
    OB = Observer(IM)
    OB2 = Observer(IM2)
    CTRL = CTRL0(IM)
    CTRL2 = CTRL0(IM2)
    sc = s_curve()
    sc2 = s_curve()
    CTRL2.pi_speed.Kp = 10
    CTRL2.pi_speed.Ti = 1.0
    # CTRL_items=vars(CTRL)
    print("Simulation Run Time={0} s".format(RUN_TIME))
    dfe = 0  # dfe for down frequency execution

    f = open(r"algorithm.dat", "w")
    fio = FileIO()
    fio.write_header_to_file(f)
    for _ in range(NUMBER_OF_LINES):

        # Command and Load Torque */
        sc.speed_ref(CTRL.timebase, 100, IM)
        sc2.speed_ref(CTRL2.timebase, 100, IM2)
        # if CTRL.timebase > 1.0:
        #     IM.Tload = 0 * IM.rpm * 0.05
        # else:
        #     IM.Tload = 0
        # Simulated ACM */
        if IM.machine_simulation(CTRL.timebase) or IM2.machine_simulation(
            CTRL2.timebase
        ):
            print("Break the loop.\n")
            break
        IM.Js += 0.5
        IM2.Js -= 0.5
        IM.Tload = 300 * (IM.x4 * IM.RAD_PER_SEC_2_RPM - IM2.x4 * IM2.RAD_PER_SEC_2_RPM)
        IM2.Tload = -300 * (
            IM.x4 * IM.RAD_PER_SEC_2_RPM - IM2.x4 * IM2.RAD_PER_SEC_2_RPM
        )
        dfe += 1
        if dfe == DOWN_FREQ_EXE:
            dfe = 0

            # Time */
            CTRL.timebase += TS
            CTRL2.timebase += TS
            OB.measurement(IM, CTRL)
            OB2.measurement(IM2, CTRL2)
            OB.observation(CTRL)
            OB2.observation(CTRL2)

            fio.write_data_to_file(f, IM, IM2)

            CTRL.control(IM.rpm_cmd, 0, OB, IM)
            CTRL2.control(IM2.rpm_cmd, 0, OB2, IM2)
        IM.inverter_model(CTRL)
        IM2.inverter_model(CTRL2)
    print("Simulation Costs Time={0} s".format(time.time() - start_time))
    f.close()

    ACMPlot2.draw_plotly(RUN_TIME, MACHINE_TS_INVERSE)


if __name__ == "__main__":
    main()
