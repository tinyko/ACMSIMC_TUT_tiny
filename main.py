from ACM import ACIM
from CTRL import CTRL0, s_curve
from Observer import Observer
from Utils import *
from FileIO import FileIO

# import Macros as mc
import ACMPlot

# import cProfile

# Debug mode raise all runtime warning
np.seterr(all="raise")


def main():
    start_time = time.time()
    IM = ACIM()
    OB = Observer(IM)
    CTRL = CTRL0(IM)
    sc = s_curve()

    CTRL_items=vars(CTRL)

    dfe = 0  # dfe for down frequency execution

    f = open(r"algorithm.dat", "w")
    fio = FileIO()
    fio.write_header_to_file(f)
    for _ in range(NUMBER_OF_LINES):

        # Command and Load Torque */
        sc.speed_ref(CTRL.timebase, IM)
        IM.Tload = 10 * np.sign(IM.rpm)  # No-load test

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
    print("Simulation time=", time.time() - start_time)
    f.close()

    ACMPlot.draw_trend()


if __name__ == "__main__":
    main()
