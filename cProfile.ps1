#Get History commands
#(Get-PSReadlineOption).HistorySavePath
python -m cProfile -o log.txt .\main.py
cprofilev -f .\log.txt