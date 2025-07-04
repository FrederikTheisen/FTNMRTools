#!/bin/tcsh

### RELAXATION DATA ANALYSIS SCRIPT ###
### Frederik Theisen, IBS, 2025
### Designed to be compatible with sparky peaklists exported from ccpnmr 3 using FT_exportSparky.py
### 	This file contains information on peak overlap that is used to flag results

### SETTINGS ###
# name of the sparky peak file.
set SPARKY_PEAKS = *.list

# The data file and format.
set DATA_FILE = ./m.ft3
set DATA_FMT = nmrpipe

# The kind of relaxation data.
set DATA_TYPE = r2 #r1, r2, noe


echo '#################'
echo "###### FIT ######"
echo '#################'

# name of the sparky peak file.
set PEAK_FMT = sparky3d

# File containing fit header and estimated uncertainties.
set RLX_DFIT_HEADER = relax_header.txt
set RLX_DFIT_ESD = vEstNoise.txt

# Names of temporary files.
set DATA_XTR_OUT = data_xtr.out
set DATA_XTR_SORT_OUT = data_xtr_sort.out
set RLX_DFIT_IN = rlx_dfit.inp

# Final output file.
set RLX_DFIT_OUT = rlx_dfit.out
set DATA_NAME = ($SPARKY_PEAKS:as/./ /)
echo $DATA_NAME
set DATA_OUT = $DATA_NAME[1]

#################################################################

echo 'Extracting data from spectrum...'
data_xtr   -data  $DATA_FILE \
     -dfmt  $DATA_FMT \
     -peak  $SPARKY_PEAKS \
     -pfmt  $PEAK_FMT \
     >  $DATA_XTR_OUT

#echo  'Sort the output.'
#sort $DATA_XTR_OUT > $DATA_XTR_SORT_OUT
echo 'Done'
echo
echo "Least-squares fit of $DATA_TYPE data..."

# Put on a fit header.
cat $RLX_DFIT_HEADER $DATA_XTR_OUT > $RLX_DFIT_IN

# Branch out according to data type.
if ($DATA_TYPE == 'r1') then
   rlx_dfit  -mode exp2 -rms_noise \
      -data $RLX_DFIT_IN \
      -d_esd $RLX_DFIT_ESD \
      -r_init  1.0 \
      -out $RLX_DFIT_OUT
else if ($DATA_TYPE == 'r2') then
   rlx_dfit  -mode exp2 -rms_noise \
      -data $RLX_DFIT_IN \
      -d_esd $RLX_DFIT_ESD \
      -r_init  10.0 \
      -out $RLX_DFIT_OUT
else if ($DATA_TYPE == 'noe') then
   rlx_dfit  -mode noe -rms_noise \
      -data $RLX_DFIT_IN \
      -d_esd $RLX_DFIT_ESD \
      -out $RLX_DFIT_OUT
endif

echo 'Done'

echo
echo 'Reorganizing output files...'

rm results.out
set infile = rlx_dfit.out

# Get the total number of lines in the file.
set numlines = `wc -l < $infile`
@ i = 1
while ( $i <= $numlines )
    set noglob
    set line = `sed -n "${i}p" $infile`
    # Split the line into fields.
    set fields = ( $line )
    unset noglob
    # If the first field is not "#", process and print it.
    if ("$fields[1]" != "#") then
        # Remove the first two characters from the first field.
        set modfield = `echo "$fields[1]" | sed "s/^..//"`
        echo "$modfield $fields[2] $fields[3]" >> pre_results.tmp

    endif
    @ i++
end

# Remove letters from assingment name
awk '{ match($1, /[A-Z]([0-9]+)[A-Z]-/, arr); print arr[1], $2, $3 }' pre_results.tmp > "${DATA_OUT}_results.out"

# Make assignment overlap map
set idx = -1
set idx = `awk 'NR == 1 {for(i=1;i<=NF;i++) { if($i=="overlap") { print i } } }' $SPARKY_PEAKS`
echo 'overlap column idx: '$idx
awk 'NR > 1 {print $1, $'$idx'}' $SPARKY_PEAKS > overlap_map.tmp

# Make second results file with overlapping peaks data
awk 'FNR==NR {keep[FNR]=$2; next} keep[FNR]=="True"' overlap_map.tmp "${DATA_OUT}_results.out" > "${DATA_OUT}_results_peak_overlap.out"
awk 'FNR==NR {keep[FNR]=$2; next} keep[FNR]=="False"' overlap_map.tmp "${DATA_OUT}_results.out" > "${DATA_OUT}_results_filtered.out"

rm plot.gnu
rm intensity.out

set file1 = relax_header.txt
set file2 = data_xtr.out

# Process file1: read the last line and extract the number of x-values and the x-values themselves
set last_line = `tail -1 $file1`
# Split the line into an array (space-separated)
set file1_fields = ( $last_line )
# The first field is the number of x-values (assumed to be an integer)
set numcols = $file1_fields[1]

# Store the x-values (fields 2 through numcols+1) in an array xvals
set xvals = ()
@ j = 1
while ($j <= $numcols)
    @ idx = $j + 1
    set xvals = ( $xvals $file1_fields[$idx] )
    @ j = $j + 1
end

# Write header commands to plot.gnu
cat << EOF > plot.gnu
#set border 31 lw 3
#set nokey
set term pdf enhanced color solid size 6,4
set output "$SPARKY_PEAKS\_intensity.pdf"
set xlabel '{Delay (s)}'
set ylabel '{Intensity (a.u.)}'
set grid
EOF

# --- Process file2 line by line
# Get the number of lines in file2
set n_res = `cat $file2 | wc -l`
set n_res = $n_res[1]
@ i = 1
while ( $i <= $n_res )
    # Extract the i-th line from file2
    set noglob
    set line = `sed -n "${i}p" $file2`
    # Split the line into fields
    set fields = ( $line )
    set label = $fields[1]
    set idx = `expr $i - 1`
    set overlap_state = `grep -w "^$label" overlap_map.tmp | awk '{print $2}'`

    if ( "$overlap_state" == "True" ) then 
    set label = "$label*" 
    endif
    
    unset noglob
    
    # Write comment lines: using the first field (label)
    echo "#$label index $idx" >> intensity.out
    echo "#$label index $idx" >> plot.gnu

    # Loop over j from 1 to numcols."
    @ j = 1
    while ( $j <= $numcols )
        # In file2 the value to plot is in field j+1
        @ idx2 = $j + 1
        set numerator = `printf "%f" $fields[$idx2]`
        # Normalize by field 2 (norm column is fixed to 2)
        set denominator = `printf "%f" $fields[2]`
        # Calculate ratio (using bc for floating-point division)
        set ratio = `echo "scale=4; $numerator / $denominator" | bc`
        set ratio = `printf "%f" $ratio`
        # Print: x-value from file1 and the computed ratio
        echo "$xvals[$j] $ratio" >> intensity.out
        @ j = $j + 1
    end
    
    # Append two blank lines to separate data blocks. These are necessary for gnuplot to recognize the new data block...
    echo "  " >> intensity.out
    echo "  " >> intensity.out
    
    # Append the fitting commands for this block
    echo "f(x)=a*exp(-b*x)" >> plot.gnu
    echo "a=1" >> plot.gnu
    echo "b=4" >> plot.gnu
    echo "fit f(x) 'intensity.out' index $idx using 1:2 via a, b" >> plot.gnu
    echo "plot [0:2.0] [-0.2:1.2] 'intensity.out' index $idx using 1:2 pt 7 t '$label', f(x) lw 3" >> plot.gnu

    @ i = $i + 1
end

# Append final plot commands
echo "set xlabel '{Residue}'" >> plot.gnu
echo "set ylabel '{$DATA_TYPE}'" >> plot.gnu
echo "plot  '${DATA_OUT}_results.out' using 1:2:3 wi er lw 2, '${DATA_OUT}_results.out' using 1:2 wi lines lw 3" >> plot.gnu
echo "plot  '${DATA_OUT}_results.out' using 1:2 wi impulses lw 3 lc rgb 'black' title 'Isolated peaks', '${DATA_OUT}_results_peak_overlap.out' using 1:2 wi impulses lw 3 lc rgb 'red' title 'Overlapping'" >> plot.gnu

echo 'Done'

/usr/bin/gnuplot plot.gnu

rm *.tmp


