#!/bin/csh
nmrPipe -verb -in test.fid                                   \
      |   nmrPipe  -fn SOL                                \
      |   nmrPipe  -fn GM -g1 0.0 -g2 8 -g3 0. -c 0.5            \
      |   nmrPipe  -fn ZF -auto                                  \
      |   nmrPipe  -fn FT                                        \
      |   nmrPipe  -fn PS -p0 60 -p1 0.0  -di                    \
      |   nmrPipe  -fn EXT -left -sw              \
      |   nmrPipe  -fn EXT -xn 6.5ppm -x1 11ppm -sw              \
      |   nmrPipe  -fn TP                                        \
#     |   nmrPipe  -fn LP -ord 16 -fb                            \
#      |   nmrPipe  -fn GM -g1 0.0 -g2 12 -g3 0. -c 1.00           \
      |  nmrPipe -fn SP -off 0.42 -pow 2 -end 0.99 -c 0.50 \
#      |   nmrPipe  -fn PS -ls -1ppm -sw                          \
      |   nmrPipe  -fn ZF -auto                                  \
      |   nmrPipe  -fn ZF -zf 1                                  \
      |   nmrPipe  -fn FT                                 \
      |   nmrPipe  -fn PS -p0 -88 -p1 180 -di                   \
      |   nmrPipe  -fn TP                                        \
      |   pipe2xyz -ov -verb -out data/test%03d.ft2   
#
# Collect data
xyz2pipe -in data/test%03d.ft2 \
   |   nmrPipe -ov -verb -out test.ft2

