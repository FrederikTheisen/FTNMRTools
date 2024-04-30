#!/bin/tcsh -f
# Bash script to extract assigned peak intensities from sparky peak list and HSQC slice
# Script reads number of slices recorded, converts the sparky peaklist to nmr pipe and runs the seriesTab command to extract intensities of each peak
# Place script in EXPNO folder

set datadir = .

# Make list of spectra
find $datadir/data -name '*.ft2' | sort -n > spectra

# Make NMR-pipe peak list from SPARKy peak list and NMR-pipe spectrum
stPeakList.pl $datadir/data/test001.ft2 peaks_sparky > peaks.dat

# Run seriesTab
seriesTab -in peaks.dat -out peaks.ser -list spectra -dx 1 -dy 1
