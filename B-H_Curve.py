import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

#POR
#75kW 740rpm cos_phi=0.86
#nominal magnetizing current=53.723
#nominal EMF=324.09
H1=[42.5,76.5,100.0,141.8,232.3]
B1=[60.0,90.0,100.0,108.5,115.9]
#Pinch Roll 1
#7.5kW 1440rpm cos_phi=0.85
#nominal magnetizing current=5.714
#nominal EMF=323.93
H2=[48.9,82.3,100.0,133.9,211.0]
B2=[59.9,90.0,100.0,111.5,123.1]
#Bridle Roll 1
#11kW 970rpm cos_phi=0.85
#nominal magnetizing current=11.167
#nominal EMF=318.37
H3=[44.4,77.5,100.0,137.9,219.0]
B3=[60.0,90.1,100.0,108.6,114.9]
#Bridle Roll 2
#15kW 970rpm cos_phi=0.86
#nominal magnetizing current=11.689
#nominal EMF=317.39
H4=[50.2,84.4,100.0,132.9,210.2]
B4=[60.0,90.0,100.0,112.1,123.7]
#Pinch Roll2
#5.5kW 1440rpm cos_phi=0.83
#nominal magnetizing current=4.427
#nominal EMF=320.90
H5=[43.8,83.6,100.0,137.4,216.0]
B5=[60.0,90.0,100.0,109.9,121.9]
#TR
#75kW 740rpm cos_phi=0.86
#nominal magnetizing current=56.236
#nominal EMF=324.32
H6=[44.2,77.7,100,141.4,232.1]
B6=[60,90.1,100,108.8,116.1]

df=pd.read_csv(r"C:\Users\tinym\Desktop\Matlab\ACMSIMC_TUT\_femm\codes3\M-19-Steel-BH-Curve-afterJMAGsmooth.BH",sep="\s+",header=None)

def rational(x, p, q): # pade appximate
    """ https://stackoverflow.com/questions/29815094/rational-function-curve-fitting-in-python
    The general rational function description.
    p is a list with the polynomial coefficients in the numerator
    q is a list with the polynomial coefficients (except the first one)
    in the denominator
    The zeroth order coefficient of the denominator polynomial is fixed at 1.
    Numpy stores coefficients in [x**2 + x + 1] order, so the fixed
    zeroth order denominator coefficent must comes last. (Edited.)
    """
    return np.polyval(p, x) / np.polyval(q + [1.0], x)

def rational8_8(x, p0, p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6, q7):
    '''88'''
    return rational(x, [p0, p1, p2, p3, p4, p5, p6, p7], [q1, q2, q3, q4, q5, q6, q7])
def rational8_7(x, p0, p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6):
    '''87'''
    return rational(x, [p0, p1, p2, p3, p4, p5, p6, p7], [q1, q2, q3, q4, q5, q6])
def rational7_7(x, p0, p1, p2, p3, p4, p5, p6, q1, q2, q3, q4, q5, q6):
    '''77'''
    return rational(x, [p0, p1, p2, p3, p4, p5, p6], [q1, q2, q3, q4, q5, q6])
def rational7_6(x, p0, p1, p2, p3, p4, p5, p6, q1, q2, q3, q4, q5):
    '''76'''
    return rational(x, [p0, p1, p2, p3, p4, p5, p6], [q1, q2, q3, q4, q5])
def rational7_5(x, p0, p1, p2, p3, p4, p5, p6, q1, q2, q3, q4):
    '''75'''
    return rational(x, [p0, p1, p2, p3, p4, p5, p6], [q1, q2, q3, q4])
def rational6_6(x, p0, p1, p2, p3, p4, p5, q1, q2, q3, q4, q5):
    '''66'''
    return rational(x, [p0, p1, p2, p3, p4, p5], [q1, q2, q3, q4, q5])
def rational6_5(x, p0, p1, p2, p3, p4, p5, q1, q2, q3, q4):
    '''65'''
    return rational(x, [p0, p1, p2, p3, p4, p5], [q1, q2, q3, q4])
def rational6_4(x, p0, p1, p2, p3, p4, p5, q1, q2, q3):
    '''64'''
    return rational(x, [p0, p1, p2, p3, p4, p5], [q1, q2, q3])
def rational5_5(x, p0, p1, p2, p3, p4, q1, q2, q3, q4):
    '''55'''
    return rational(x, [p0, p1, p2, p3, p4], [q1, q2, q3, q4])
def rational5_3(x, p0, p1, p2, p3, p4, q1, q2):
    '''53'''
    return rational(x, [p0, p1, p2, p3, p4], [q1, q2])
def rational4_4(x, p0, p1, p2, p3, q1, q2, q3):
    '''44'''
    return rational(x, [p0, p1, p2, p3], [q1, q2, q3])
def rational4_3(x, p0, p1, p2, p3, q1, q2):
    '''43'''
    return rational(x, [p0, p1, p2, p3], [q1, q2])
def rational3_3(x, p0, p1, p2, q1, q2):
    '''33'''
    return rational(x, [p0, p1, p2], [q1, q2])

popt, pcov =curve_fit(rational8_8,df[0],df[1])

plt.figure(figsize=(10,15))
plt.rc('font', family='SimHei', size=11, weight='bold')
plt.subplot(3,1,1)
plt.title("无取向硅钢\n 磁化曲线\n 2020.02.06 LPD")
plt.ylabel("磁感应强度 $B\ [T]$")
plt.plot(df[0][0:146],df[1][0:146],'black')
plt.plot(df[0][0:146], rational8_8(df[0][0:146], *popt),'red')
plt.legend(["实际测量描点","拟合曲线"],loc="lower right")
plt.xscale('log',basex=10)
plt.grid(b=True, which='major', color=[171/255,209/255,57/255], linestyle='-')
plt.grid(b=True, which='minor', color=[171/255,209/255,57/255], linestyle='-')

plt.subplot(3,1,2)
plt.xlabel("磁场强度 $H\ [A/m]$")
plt.ylabel("磁感应强度 $B\ [T]$")
plt.plot(df[0][0:113],df[1][0:113],'black')
plt.plot(df[0][0:113], rational8_8(df[0][0:113], *popt),'red')
plt.legend(["实际测量描点","拟合曲线"],loc="lower right")
plt.grid(color=[171/255,209/255,57/255])
#plt.xscale('log',basex=10)
#plt.grid(b=True, which='major', color=[171/255,209/255,57/255], linestyle='-')
#plt.grid(b=True, which='minor', color=[171/255,209/255,57/255], linestyle='-')

plt.subplot(3,1,3)
plt.rc('font', family='SimHei', size=11, weight='bold')
plt.plot(H1,B1,H2,B2,H3,B3,H4,B4,H5,B5,H6,B6)
plt.xlabel("励磁电流 [%]\n (100%=额定励磁电流)")
plt.ylabel("磁通量 [%]\n (100%=额定磁通=额定电动势EMF)")
plt.legend(["开卷机75kW 740rpm",
            "1号夹送辊7.5kW 1440rpm",
            "1号张力辊11kW 970rpm",
            "2号张力辊15kW 970rpm",
            "2号夹送辊5.5kW 1440rpm",
            "卷取机75kW 740rpm"],loc="lower right")
plt.grid(color=[171/255,209/255,57/255])

#plt.savefig("B-H_Curve_NonOriented.pdf")

plt.show()