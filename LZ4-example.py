import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lz4,lz4.frame
import sys
import os
from tabulate import tabulate as tbl
import timeit
import json
import zstandard as zstd

lz4lver=lz4.library_version_string()
lz4ver=lz4.VERSION
df=pd.read_csv(r"trend.dat")
nd=df.to_numpy()
sp=nd.shape
tb=nd.tobytes()

with open('test2.dat','wb') as fh:
    cctx = zstd.ZstdCompressor(level=1,threads=4)
    start2 = timeit.default_timer()
    for i in range(10**3):
        compressed2 = cctx.compress(tb)
    stop2 = timeit.default_timer()
    time2=(stop2 - start2)/10**3
    fh.write(compressed2)
with open('test2.dat','rb') as fh:
    cctx = zstd.ZstdDecompressor()
    tmp = fh.read()
    start3 = timeit.default_timer()
    for i in range(10**3):
        decompressed2 = cctx.decompress(tmp)
    stop3 = timeit.default_timer()
    time3=(stop3 - start3)/10**3
zz=zstd.ZSTD_VERSION
zzs=""
for i in range(len(zz)):
    if i == 0 :
        zzs=str(zz[i])
        continue
    zzs=zzs+"."+str(zz[i])

start = timeit.default_timer()
for i in range(10**3):
    compressed = lz4.frame.compress(tb)
stop = timeit.default_timer()
time=(stop - start)/10**3

with lz4.frame.open('test.dat', mode='wb') as fp:
    writtensize=fp.write(compressed)
with lz4.frame.open('test.dat', mode='r') as fp:
    output_data = fp.read()

start4 = timeit.default_timer()
for i in range(10**3):
    dcp=lz4.frame.decompress(output_data)
stop4 = timeit.default_timer()
time4=(stop4 - start4)/10**3

tb_origin=np.frombuffer(dcp,dtype=nd.dtype)
result=np.array_equal(tb_origin.reshape(60000,12),nd)
result2=os.path.getsize('test.dat')
result3=os.path.getsize('trend.dat')
result4=os.path.getsize('test2.dat')

dic3=lz4.frame.get_frame_info(compressed)
bs=dic3['block_size']
ty=nd.dtype.name

dic=[['LZ4 library version',lz4lver],
     ['LZ4 python module version',lz4ver],
     ['LZ4 block size',"{0}kB".format(int(np.round(bs/1024)))],
     ['LZ4 compression level',0],
     ['Zstd library version',zzs],
     ['Zstd Compression level',1],
     ['Zstd Compression Threads',4],
     ['Type of element',ty],
     ['Shape of data',sp],
     ['Original Size',"{0}kB".format(int(np.round(result3/1024)))],
     ['Compressed Size(LZ4)',"{0}kB".format(int(np.round(result2/1024)))],
     ['Compressed Size(Zstd)',"{0}kB".format(int(np.round(result4/1024)))],
     ['Compression Ratio(LZ4)',"{0:.6f}".format(result3/result2)],
     ['Compression Ratio(Zstd)',"{0:.6f}".format(result3/result4)],
     ['Compression Time(LZ4)',"{0:.6f}s".format(time)],
     ['Compression Time(Zstd)',"{0:.6f}s".format(time2)],
     ['Decompression Time(LZ4)',"{0:.6f}s".format(time4)],
     ['Decompression Time(Zstd)',"{0:.6f}s".format(time3)],
     ['Original == Compressed ?',result]]
with open('data.json', 'w') as fp:
    json.dump(dic, fp)

with open('data.json', 'r') as fp:
    dic2 = json.load(fp)
print(tbl(dic2,
          tablefmt="fancy_grid",
          headers=['Name','Value']))