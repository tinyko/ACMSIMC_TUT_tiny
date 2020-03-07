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
    OB = Observer(IM)
    CTRL = CTRL0(IM)
    sc = s_curve()
    sc.ACC = 2
    sc.DCC = 2

    # CTRL_items=vars(CTRL)
    print("Simulation Run Time={0} s".format(RUN_TIME))
    dfe = 0  # dfe for down frequency execution

    f = open(r"algorithm.dat", "w")
    fio = FileIO()
    fio.write_header_to_file(f)
    for _ in range(NUMBER_OF_LINES):

        # Command and Load Torque */
        sc.speed_ref(CTRL.timebase, 100, IM)
        if CTRL.timebase > 1.0:
            IM.Tload = 0 * IM.rpm * 0.05
        else:
            IM.Tload = 0
        # Simulated ACM */
        if IM.machine_simulation(CTRL.timebase):
            print("Break the loop.\n")
            break
        dfe += 1
        if dfe == DOWN_FREQ_EXE:
            dfe = 0

            # Time */
            CTRL.timebase += TS

            OB.measurement(IM, CTRL)

            OB.observation(CTRL)

            fio.write_data_to_file(f, IM, CTRL)

            CTRL.control(IM.rpm_cmd, 0, OB, IM)

        IM.inverter_model(CTRL)
    print("Simulation Costs Time={0} s".format(time.time() - start_time))
    f.close()

    ACMPlot2.draw_plotly(RUN_TIME, MACHINE_TS_INVERSE)


if __name__ == "__main__":
    main()
