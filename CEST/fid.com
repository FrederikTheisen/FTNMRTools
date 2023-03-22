#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1600 -dspfvs 20 -grpdly 68  \
  -xN              1280  -yN                70  -zN               128  \
  -xT               614  -yT                35  -zT                64  \
  -xMODE            DQD  -yMODE           Real  -zMODE  Echo-AntiEcho  \
  -xSW        12500.000  -ySW         6034.929  -zSW         2189.323  \
  -xOBS         800.134  -yOBS         201.205  -zOBS          81.086  \
  -xCAR          4.7589  -yCAR         58.0747  -zCAR        119.5967  \
  -xLAB              HN  -yLAB             13C  -zLAB             15N  \
  -ndim               3  -aq2D         States                          \
 | nmrPipe  -ov -verb -out test.fid

nmrPipe -in test.fid \
 | nmrPipe -fn TP -auto \
 | nmrPipe -fn ZTP  \
 | nmrPipe -fn TP -auto \
 | nmrPipe -ov -verb -out tt.fid
 mv tt.fid test.fid
