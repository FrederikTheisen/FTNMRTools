#!/bin/tcsh

foreach i (data/*.ft2)
sethdr $i -ndim 2 -zN 1 -zT 1 -zMODE Real
end
