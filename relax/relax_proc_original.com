#!/bin/csh

#T1 block length in seconds
set delay_block = 0.01
set delay_constant = 0.01
set path2vc = "vclist"
set path2vd = "vdlist"
set path2ncyc = "ncyc_list"

#Modify calculation to fit pulse sequence
if (-e "$path2vc") then
   cat "$path2vc" | awk 'NF {print '${delay_constant}'+$1*'${delay_block}'}' > delaylist
else if (-e "$path2vd") then
   cat "$path2vd" | awk 'NF {unit=($1 ~ /us$/) ? 1000000 : ($1 ~ /m$/) ? 1000 : 1; print substr($1, 1, length($1)-1)/unit}' > delaylist
else if (-e "$path2ncyc") then
   cat "$path2ncyc" | awk 'NF {print $1}' > delaylist
endif

echo '#########################'
echo "### MAKE RELAX HEADER ###"
echo '#########################'

awk '{a[NR]=$0} END {printf "%d", NR; for (i=1; i<=NR; i++) printf " %s", a[i]; print ""}' delaylist > relax_header.txt

echo
echo '##########################'
echo "### RUNNING FID SCRIPT ###"
echo '##########################'

# This script can be generated using the bruker command
# Paste in the script here and do the phasing below in the proc part

bruk2pipe -verb -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 2760 -dspfvs 20 -grpdly 67.9857330322266  \
  -xN              1280  -yN              2200  \
  -xT               640  -yT              1100  \
  -xMODE            DQD  -yMODE  Complex  \
  -xSW         7246.377  -ySW         1515.152  \
  -xOBS         600.483  -yOBS          60.853  \
  -xCAR           4.771  -yCAR         118.577  \
  -xLAB              HN  -yLAB             15N  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -ov

echo
echo '#####################################' 
echo "### RUNNING SPECTRUM SPLIT SCRIPT ###"
echo '#####################################'

if (-e test.fid) then
   /bin/rm -rf fid

   mkdir fid

   set vList = (`cat delaylist`)
   set n     = $#vList
   set i     = 1
   set cList = ""

   foreach v ($vList)
      set cList = ($cList 0.0)
   end

   foreach v ($vList)
      set cList[$i] = 1.0
      set tau       = (`echo $v | sed 's/m//'`) 
      set outName   = (`printf fid/test%03d.fid $i`)

      echo $tau $outName $cList

      nmrPipe -in test.fid \
      | nmrPipe -fn COADD -time -axis Y -cList $cList -verb      \
      | nmrPipe -fn MAC -macro $NMRTXT/bruk_ranceY.M -noRd -noWr \
      -out $outName -ov

      set cList[$i] = 0.0

      @ i++
   end

   /bin/rm -rf test.fid 

   series.com fid/test*.fid

endif

echo
echo '###########################'
echo "### RUNNING PROC SCRIPT ###"
echo '###########################'

set out = proc.ucsf
set tauValues = (`cat delaylist | sed 's/m//'`)
set idxExpmt  = ( `seq 1 ${#tauValues}` )

rm -rf ft
mkdir ft

foreach n ( `seq 1 ${#tauValues}` )
  echo Time t${n} \($tauValues[$n]\)...
  set outName = (`printf ft/test%03d.ft2 $idxExpmt[$n]`)
  set inName  = (`printf fid/test%03d.fid $n`)
  echo "* input  = $inName"
  echo "* output = $outName"

  nmrPipe -in $inName \
  | nmrPipe  -fn SOL                             \
  | nmrPipe  -fn SP -off 0.4 -end 0.98 -pow 2 -c .5    \
  | nmrPipe  -fn ZF -size 4096                          \
  | nmrPipe  -fn FT                                     \
#  | nmrPipe  -fn EXT -x1 11ppm -xn 6ppm -sw -verb       \
  | nmrPipe  -fn EXT -left -sw                               \
  | nmrPipe  -fn PS -p0 0 -p1 0 -di                     \
  | nmrPipe  -fn TP                                     \
  | nmrPipe  -fn SP -off 0.4 -end 0.98 -pow 2 -c .5    \
  | nmrPipe  -fn ZF -size 1024                           \
  | nmrPipe  -fn FT                                      \
  | nmrPipe  -fn PS -p0 0 -p1 0 -di                 \
  | nmrPipe  -fn TP                                     \
  | nmrPipe  -fn POLY -auto -ord 4                      \
     -verb -ov -out $outName

  #@ i = $n - 1
  sethdr $outName -tau $tauValues[$n]
end

series.com ft/test*.ft2

xyz2pipe -in ft/test%03d.ft2 > m.ft3

pipe2ucsf m.ft3 $out
ucsfdata -a1 15N -o1 10. -sw1 100. -f1 100. -a2 15N $out

# Convert 2D planes to 2D spetra for ccpnmr analysis viewing/peak picking

foreach i (ft/*.ft2)
  sethdr $i -ndim 2 -zN 1 -zT 1 -zMODE Real
  setfdata $i -fileCount 1
end

echo
echo '######################'
echo "### Estimate Noise ###"
echo '######################'

# Initialize variables
set noiseList = ()
set count = 0

# Loop over files in ./ft directory
foreach spec (ft/test*.ft2)
  echo $spec
  set output = `specStat.com -in $spec -stat vEstNoise`
  echo $output
  set noise = $output[2]
  set noiseList = ($noiseList $noise)
  @ count++
end

echo "$count $noiseList" > vEstNoise.txt



