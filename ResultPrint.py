from tabulate import tabulate as tbl
import json
import sys

import subprocess

p = subprocess.Popen(
    "powershell.exe  Get-ComputerInfo -Property @('CsProcessors','CsSystemFamily','WindowsBuildLabEx','OsName','OsVersion') ",
    stdout=subprocess.PIPE,
    shell=True,
)
(output, err) = p.communicate()

print(output.decode("gb2312"))
path = sys.argv[1]
with open(path, "r") as fp:
    dic2 = json.load(fp)
print(tbl(dic2, tablefmt="fancy_grid", headers=["Name", "Value"]))
