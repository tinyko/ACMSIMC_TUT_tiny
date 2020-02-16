## Tiny 2020.02.13(2020.02.16 updated)
TO DO List
1) ~~Seperate Programm Functions into different python modules~~
2) GUI in QT (Auto Scale GUI, QT GUI layout)
3) Async communications between Real-time and GUI module
4) design a python to c source translator(like simulink) (further plan)
5) motor topology
6) add json config file to manage parameters
7) create a **Trend Viewer** to evaluate trends data , which could have some functions to calculate,zoom in/out,.etc(like IBA©analyzer)
8) save all **critical** data in trends
9) create ~~a speical~~~ **HDF5** data file structure to contain data informations(timebase,timestamp,.etc) and values(bool,string,int,float,.etc)

Done List
1) Converted C program into Python Single File "Alter"
2) Modified ACMPlot as a library to Main Python Script
3) Compared LZ4 and Zstandard (Winner is LZ4! In Realtime scenanio)
4) numpy ndarray could be converted into bytes , vice versa
