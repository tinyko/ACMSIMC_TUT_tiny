from Utils import *


class FileIO:
    def __init__(self):
        self.j = 0

    def write_header_to_file(self, f):
        f.write(
            "ACM.x[0],ACM.x[1],ACM.x[2],ACM.x[3],ACM.x[4],ACM.Tem,CTRL.uMs_CMD,CTRL.uTs_CMD,CTRL.iMs_e,CTRL.iTs_e,ACM.rpm_cmd,e_omega\n"
        )

        f2 = open(r"info.dat", "w")
        f2.write("TS,DOWN_SAMPLE\n")
        f2.write("{0}, {1}\n".format(TS, DOWN_SAMPLE))
        f2.close()

    def write_data_to_file(self, f, IM, CTRL):
        self.j += 1
        if self.j == DOWN_SAMPLE:
            self.j = 0
            f.write(
                "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n".format(
                    IM.x[0],
                    IM.x[1],
                    IM.x[2],
                    IM.x[3],
                    IM.x[4],
                    IM.Tem,
                    CTRL.uMs_cmd,
                    CTRL.uTs_cmd,
                    CTRL.iMs - CTRL.iMs_cmd,
                    CTRL.iTs - CTRL.iTs_cmd,
                    IM.rpm_cmd,
                    IM.rpm_cmd - IM.x[4] * IM.RAD_PER_SEC_2_RPM,
                )
            )
