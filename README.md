# FTNMRTools

Various tools to handle NMR processing and analysis tasks.
To download a file, open the file and right click the 'Raw' button in top right of the file area and click download archive.

## CEST

Collection of tools and files needed to process and analyze CEST data in SBiNLab.

## 1DBaselineCorrection:

Python tool to baseline correct data assuming a polynomial baseline. The script should work on MacOS, Linux and Windows.

Will also let you integrate regions of the loaded spectra.

The script can read two file types:
1) Topspin exported file in the following format:

```
# File created = Wednesday, March 8, 2023 2:44:20 PM CET
# Data set = expname  121  1  /path/to/folder
# Spectral Region:
# F1LEFT = 0.0 s. F1RIGHT = 0.0049992 s.
# F2LEFT = 2.8 ppm. F2RIGHT = 0.0 ppm.
#
# NROWS = 16 ( = number of points along the F1 axis)
# NCOLS = 5874 ( = number of points along the F2 axis)
...
# row = 0
1.55127808E9
1.536176128E9
1.522184192E9
1.525284864E9
...
```

2) Topspin 'totxt' exported files in a folder. 'multitotxt' can be executes by downloading the AU program from here and running that in Topspin

Topspin command line
```
multitotxt
```
The AU program will ask for the number of experiments and the output path. The output files will be named according to the EXPNO of the experiments.

Usage of the script:
```
#If all data sets are in the same file
python3 1DBaselineCorrection.py path-to-file
#If multiple files exported using multitotxt (or totxt many times)
python3 1DBaselineCorrection.py path-to-folder-with-files
```
Output is 2 files containing the baseline corrected data and the integrated peaks.

baselined.txt: [offset in column 1 with corrected data in subsequent columns]
```
ppm 0 1 ...
2.8 8871249.127078533 38488362.885465145 ...
```

peaks.txt: [peak# in column 1 and subsequent columns contain the areas of the peaks for each spectra]
```
peak 0 1 ... 
0 97437044817.22278 89422449114.04411 ...
```
## GetOffsetTimes:

Python tool to get experiment times from a folder containing experiments. Script is tested on Linux.

Tool is designed to collect offset times for series of experiments where the first experiment of each series is number x001. The script will look for folders in the currect directory.
Folder structure should be as follows:
```
1         # exp is ignored
2         # exp is ignored
1001      # series 1 reference exp and first data point
1002      # s1 point 2
...
2001      # series 2 reference exp and first data point
2002      # s2 point 2
...
```

Usage:
```
# navigate to folder with spectra
# save GetOffsetTimes.py to this folder

python3 GetOffsetTimes.py
```

Output: Console dump
```
Exp #1001, time: 1678887925
1001 0
1002 65
1003 130
...
```

------
