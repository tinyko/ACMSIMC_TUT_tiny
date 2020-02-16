# AC Induction Motor Simulation in Python
![Screenshot](./Screenshot/SeperateModules.jpg "Screenshot in vscode")

## Introduction
This project aims to simulate <a href="https://en.wikipedia.org/wiki/Induction_motor">AC induction motor</a> in Python.

* Forked from https://github.com/horychen/ACMSIMC_TUT

* Thanks to <a href="https://horychen.github.io/about/">Jiahao Chen | 陈嘉豪</a>

* GPL-3.0 licensed

* Code structure
   * **Main Loop**
      * **Main** : Simulation Main Loop
      * **ACM** : AC Induction Motor Class
      * **CTRL** : Controller Class
      * **Observer** : Observer Class
      * **Utils** : All Utils imports and constants
      * **FileIO** : File Input and Output class
      * **Macros** : All lambda functions
      * **ACMPlot** : Draw Treand plot
   * **Addition Files**
      * **B-H_Cruve** : Draw B-H curve from 33.txt
      * **LZ4-example** : Test LZ4 and Zstandard performance
      * **ResultPrint** : Print result table in tabulate

* Prerequisites
   * Python V3.8
   * 
## Tiny 2020.02.16
TO DO List
1) ~~Seperate Programm Functions into different python modules~~
2) GUI in PyQT5 (~~Auto Scale GUI, QT GUI layout~~)
3) Async communications between Real-time and GUI module
4) Design a python to c source translator(like simulink) (further plan)
5) ~~motor topology~~
6) Add json config file to manage parameters
7) create a **Trend Viewer** to evaluate trends data , which could have some functions to calculate,zoom in/out,.etc(like IBA©analyzer)
8) Save all **critical** data in trends
9) Create ~~a speical~~ **HDF5** data file structure to contain data informations(timebase,timestamp,.etc) and values(bool,string,int,float,.etc)
10) Remove all **redundant** code and cross-reference
11) Seperate ODE-solver to a single module
12) Decoupling sub-modules

Done List
1) ~~Converted C program into Python Single File "Alter"~~
2) ~~Modified ACMPlot as a library to Main Python Script~~
3) ~~Compared LZ4 and Zstandard (Winner is LZ4! In Realtime scenanio)~~
4) ~~numpy ndarray could be converted into bytes , vice versa~~
6) **addition files**:
   * **B-H_Cruve** : Draw B-H curve from 33.txt
   * **LZ4-example** : Test LZ4 and Zstandard performance
   * **ResultPrint** : Print result table in tabulate
