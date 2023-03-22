#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1600 -dspfvs 21 -grpdly 76  \
  -xN              2048  -zN               200  -yN                84  \
  -xT              1024  -zT               100  -yT                84  \
  -xMODE            DQD  -zMODE  Echo-AntiEcho  -yMODE           Real  \
  -xSW        12500.000  -zSW         2189.321  -ySW          810.860  \
  -xOBS         800.134  -zOBS          81.086  -yOBS          81.086  \
  -xCAR           4.7750  -zCAR         119.5947  -yCAR         117.078  \
  -xLAB              HN  -zLAB            15Ny  -yLAB            15Nz  \
  -ndim               3  -aq2D         States                         \
 | nmrPipe  -ov -verb -out test.fid

nmrPipe -in test.fid \
 | nmrPipe -fn TP -auto \
 | nmrPipe -fn ZTP  \
 | nmrPipe -fn TP -auto \
 | nmrPipe -ov -verb -out tt.fid
 mv tt.fid test.fid
