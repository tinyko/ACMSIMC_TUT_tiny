from Utils import *


class FileIO:
    def __init__(self):
        self.j = 0

    def write_header_to_file(self, f):
        f.write(
            "ACM.x0,ACM.x1,ACM.x2,ACM.x3,ACM.x4,ACM.Tem,ACM.ual,ACM.ube,ACM.ial,ACM.ibe,ACM.rpm_cmd,e_omega\n"
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
                    IM.x0,
                    IM.x1,
                    IM.x2,
                    IM.x3,
                    IM.x4,
                    IM.Tem,
                    IM.ual,
                    IM.ube,
                    IM.ial,
                    IM.ibe,
                    IM.rpm_cmd,
                    IM.rpm_cmd - IM.x4 * IM.RAD_PER_SEC_2_RPM,
                )
            )
