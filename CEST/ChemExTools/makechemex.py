import sys
import re
import os
from os import listdir
from os.path import isfile, join

### Frederik Theisen 2022 ###
### Generates ChemEx project folder structure based on processed CEST data
### To run this script, you should prepare you files according to the instructions
### You should also fill out the SETTINGS part of the script
### Make an EXPERIMENT folder named according to the experiment: eg 15N
### Inside the EXPERIMENT folder, make and INPUT folder which contains the intensity files and the frqlist file
###### peaks.ser files should be name B1.ser, where B1 is the B1 field and should be written as eg 12p5.ser or 12_5.ser
### Place the makechemex.py script in the parent folder of the EXPERIMENT folder.
###### Example Folder Structure:
# /EXPERIMENT1/
# //INPUT/
# ///12p5.ser
# ///6p25.ser
# ///25.ser
# ///frq1list
# /EXPERIMENT2/
# //INPUT/
# ///25.ser
# ///frq3list
# /makechemex.py
### Navigate a terminal to the EXPERIMENT# folder and run the makechemex.py script using the INPUT folder as argument like this: python3 ../makechemex.py ./INPUT

### HOW TO RUN SCRIPT ###
### In the terminal write the following: python3 makechemex.py <path-to-INPUT-folder>
### You should replace <path-to-INPUT-folder> with the path to the INPUT folder
### You can do this by writing "python3 makechemex.py " in the terminal and then dragging the INPUT folder onto the terminal a letting go, this will paste the path to the dragged folder
### Click 'return'/'enter' to execute the line (the button that makes a new line in word)

### Instructions for preparing makechemex.py script
### You should only change values in the SETTINGS part of the script
### Get B0 from reference experiment
### Get OFFSET from B1 dimension offset in acqupars
### MEASURED_OFFSET is obtained by running calcOffset in the CEST experiment

#SETTINGS
EXPERIMENTTYPE = "cest_15n"	#Experiment type. Check ChemEx github for options. Defaults are 'cest_15n' and 'cest_13c'.
PB = 0.042			#Initial guess percent bound
KEX = 350			#Initial guess Kex
TAU = 10 			#R2 initial guess related parameter
T1TIME = 0.4 			#CEST delay time (for 15N CEST probably = d21)
OFFSET = 119.5			#B1 reference ppm used to obtain B1 offset Hz 
MEASURED_OFFSET = 119.12345	#Referenced B1 offset ppm used output correct chemical shifts
B0 = 800.134			#Proton B0 field strength MHz
DLABEL13C = True		#Is sample also 13C labelled

#AUTOASSIGNED SETTINGS 		#Set to 'None' for auto
INPUTTYPE = None 		#ppm or hz
CESTATOM = None 		#CEST atom label
CESTATOMB0 = None		#CEST atom larmor frq MHz, used to convert ppm B1 offsets to Hz. Only used if input is in ppm


#OUTPUT DIRECTORY STRUCTURE
PATH = sys.argv[1]
OUTPUTFOLDER = "PROJECT"
DATAFOLDER = "DATA"
PROJECTFOLDER = "CONFIGURATION"
EXPERIMENTFOLDER = "EXPERIMENTS"
PARAMETERFOLDER = "PARAMETERS"
METHODFOLDER = "METHOD"
RESULTFOLDER = "RESULTS"

def getfile(filename):
	return PATH + "/" + filename

def getoutputpath(folder):
	return "./" + OUTPUTFOLDER + "/" + folder

def getprojectfolder(folder):
	return getoutputpath(PROJECTFOLDER) + "/" + folder

def satfrqoffset(ppm):
	dppm = ppm - OFFSET
	return dppm

def ppmtohz(ppm):
	return ppm * CESTATOMB0

def writeline(file, line):
	file.write(line + "\n")

cestatoms = {'cest_15n': 'N', 'cest_13c':'C'}
files = [f for f in listdir(PATH) if isfile(join(PATH,f))]
serfiles = [f for f in files if ".ser" in f]
listfile = [f for f in files if "list" in f][0]
badassigns = []

if CESTATOM is None: CESTATOM = cestatoms[EXPERIMENTTYPE]
if CESTATOMB0 is None:
	if CESTATOM == 'N': CESTATOMB0 = B0 / 9.8650113766
	elif CESTATOM == 'C': CESTATOMB0 = B0 / 3.9760539474

print("Data files:   ",serfiles)
print("List file:    ",getfile(listfile))
print("Spectrometer: ",B0," T1 = ", str(T1TIME)+"s")
print("CEST atom:    ",CESTATOM,CESTATOMB0)
print("B1 OFFSET:    ",OFFSET)

satfrqs = []
satfrqcount = 0

inputtype = INPUTTYPE #set input type (if 'None', then program will try to guess)

with open(getfile(listfile)) as lf:
	lines = lf.readlines()

	if INPUTTYPE == None: #set file specific input type
		if lines[0].strip() == "bf ppm": inputtype = 'ppm'
		else: inputtype = 'hz'

	print("FORMAT: " + str(inputtype))

	for line in lines[1:]:
		line = line.strip()
		if line == "": continue
		elif line == "bf ppm": continue
		frq = float(line)
		if inputtype == 'ppm': satfrqs.append(ppmtohz(satfrqoffset(frq)))
		elif inputtype == 'hz': satfrqs.append(frq)
		satfrqcount += 1
###
### PREPARE OUTPUT DATA STRUCTURE
###
if not os.path.exists(getoutputpath("")): os.mkdir(getoutputpath(""))
if not os.path.exists(getoutputpath(DATAFOLDER)): os.mkdir(getoutputpath(DATAFOLDER))
if not os.path.exists(getoutputpath(PROJECTFOLDER)): os.mkdir(getoutputpath(PROJECTFOLDER))
if not os.path.exists(getoutputpath(RESULTFOLDER)): os.mkdir(getoutputpath(RESULTFOLDER))
if not os.path.exists(getprojectfolder(EXPERIMENTFOLDER)): os.mkdir(getprojectfolder(EXPERIMENTFOLDER))
if not os.path.exists(getprojectfolder(PARAMETERFOLDER)): os.mkdir(getprojectfolder(PARAMETERFOLDER))
if not os.path.exists(getprojectfolder(METHODFOLDER)): os.mkdir(getprojectfolder(METHODFOLDER))

###
### PRODUCE RESIDUE DATA FILES FROM .SER FILES
###

allresidues = []
experiments = {}

for file in serfiles:
	B1 = float(file.split(".")[0].replace('p','.').replace('_','.'))
	expdir = getoutputpath(DATAFOLDER + "/" + str(B1).replace(".","_"))
	if not os.path.exists(expdir): os.mkdir(expdir)
	with open(getfile(file)) as f:
		lines = f.readlines()
		badass_idx = 0

		residues = {}

		for line in lines:
			dat = line.strip().split()
			if len(dat) < 2: continue
			if dat[0].isdigit():
				residue = {}
				residue["ppm"] = float(dat[4])
				ass = re.sub('[^0-9]','',dat[6].split('-')[0])
				residue_type = dat[6][0]
				if ass == "": #Handle unassigned residues
					ass = dat[3].replace('.','') #Residue number set to N ppm
					#if len(badassigns) <= badass_idx:
					#	ass = input("Assign residue unique number ["+str(badass_idx)+"] (" + dat[6] + "): ")
					#	badassigns.append(ass) 								#Store in array for auto assign in subsequent files
					#else: ass = badassigns[badass_idx]
					#badass_idx += 1
				resi = int(ass)
				residue["filename"] = dat[6].split('-')[1]
				residue["i"] = dat[8:]
				while resi in residues:
					print(resi)
					resi = resi + 1000
				tmpname = residue["filename"][1:]
				while tmpname[0].isdigit(): tmpname = tmpname[1:]
				residue["residue_name_id"] = residue_type + str(resi) + tmpname
				residues[resi] = residue

		for resi in residues:												#produce residue files
			residue = residues[resi]
			if resi not in allresidues: allresidues.append(resi)
			with open(expdir + "/" + residue["filename"] + ".out", "w") as f:
				f.write("#Offset (Hz)\tIntensity\tUncertainty\n")
				writeline(f,"  -1.000e+05   1  0.1")									#Necessary reference value.
				for i in range(satfrqcount):
					f.write("  " + str(satfrqs[i]) + "\t" + str(residue["i"][i]) + "\t" + str(0.1) + "\n")	#data
	experiments[B1] = {"expdir": str(B1).replace(".","_"), "res":residues}


#Order the list of all residues
allresidues = sorted(allresidues)

###
### MAKE .TOML CHEMEX PROJECT FILES
### ONE FILE IS CREATED FOR EACH B1
### MULTIPLE FILES ARE LOADED USING: -e file1 file2 ...
###
for b1 in experiments:
	exp = experiments[b1]

	with open(getprojectfolder(EXPERIMENTFOLDER) + "/" + str(b1) + ".toml", "w") as f:
		writeline(f, "[experiment]")
		writeline(f, "name         = \"" + EXPERIMENTTYPE + "\"")
		writeline(f, "time_t1      = " + str(T1TIME))
		writeline(f, "carrier      = " + str(MEASURED_OFFSET)) #Replaced OFFSET with MEASURED_OFFSET
		writeline(f, "b1_frq       = " + str(b1))
		writeline(f, "b1_inh_scale = inf") #B1 inhomogeneity expressed as a fraction of 'b1_frq'.
		writeline(f, "")
		writeline(f, "[conditions]")
		writeline(f, "h_larmor_frq = " + str(B0))
		if DLABEL13C: writeline(f, "label = [\"13C\"]")
		writeline(f, "")
		writeline(f, "[data]")
		writeline(f, "path           = \"../../DATA/" + exp["expdir"] + "/\"")
		writeline(f, "error          = \"scatter\"")
		#writeline(f, "filter_offsets = [[0.0, 26.0]] #unknown importance and not sure how this filter is implemented (in ppm)
		writeline(f, "")
		writeline(f, "  [data.profiles]")
		for resi in exp["res"]:
			res = exp["res"][resi]
			writeline(f, "  " + res["residue_name_id"] + " = \"" + res["filename"] + ".out\"")

###
### MAKE .TOML CHEMEX PARAMETERS FILE
### PROVIDES HSQC 15N OFFSETS FOR INITIAL GUESSING OR FIXED VALUE
### DOES NOT PROVIDE INITIAL GUESSES FOR MINOR POPULATION OFFSETS
### USE chemex pick_cest -e file -o outputfolder FOR THIS
###

parameterfile = getprojectfolder(PARAMETERFOLDER) + "/parameters.toml"
if not os.path.exists(parameterfile):
	with open(parameterfile, "w") as f:
		writeline(f, "[GLOBAL]")
		writeline(f, "PB     = " + str(PB))
		writeline(f, "KEX_AB = " + str(KEX))
		writeline(f, "TAUC_A = " + str(TAU))
		if CESTATOM == 'N':
			writeline(f, "")
			writeline(f, "[CS_A]")
			residuelist = {}
			for b1 in experiments:
				residuelist = experiments[b1]["res"]
				break
			for resi in allresidues:
				res = residuelist[resi]
				writeline(f, res["residue_name_id"] + " = " + str(res["ppm"]))
else: print("PARAMETER file already exists.")

###
### MAKE .TOML METHOD FILE
### FILE DEFINES WHICH PARAMETERS ARE FITTED AND WHICH ARE FIXED
###

methodfile = getprojectfolder(METHODFOLDER) + "/method.toml"
if not os.path.exists(methodfile):
	with open(methodfile, "w") as f:
		writeline(f, "[METHOD]")
		writeline(f, "FIT = [\"CS_A\"]") #Do not fix the HSQC specified peak position
else: print("METHOD file already exists.")

methodfile = getprojectfolder(METHODFOLDER) + "/globalmethod.toml"
if not os.path.exists(methodfile):
	with open(methodfile, "w") as f:
		writeline(f,"""[STEP1]
INCLUDE = [1,2,3]		
FIT     = ["CS_A"]

[STEP2]
INCLUDE = "ALL"
FIX     = ["PB", "KEX_AB"]

[STEP3]
FIT     = ["PB", "KEX_AB"]""")
else: print("GLOBALMETHOD file already exists.")

with open("./" + OUTPUTFOLDER + "/guide.sh","w") as f:
	writeline(f,"#!/bin/bash")
	writeline(f,"")
	writeline(f,"#Open terminal in parent directory. Parent dir contains the 'PROJECT' folder.")
	writeline(f,"#Activate Chemex conda environment")
	writeline(f,"conda activate chemex")
	writeline(f,"")
	writeline(f,"#Command to view data and pick minor states")
	writeline(f,"chemex pick_cest -e " + getprojectfolder(EXPERIMENTFOLDER) + "/*.toml -o " + getprojectfolder(PARAMETERFOLDER))
	writeline(f,"#Go to the new sandbox directory in the parent dir and copy the [DW_AB] values to your parameters file.")
	writeline(f,"#Rerun this program if you need to pick multiple states. Then manually rename the [DW_AX] variable name.")
	writeline(f,"")
	writeline(f,"#Fit all profiles using a no exchange model (used to evaluate if profiles have exchange or not)")
	writeline(f,"chemex fit -e " + getprojectfolder(EXPERIMENTFOLDER) + "/*.toml"
		+ " -p " + getprojectfolder(PARAMETERFOLDER) + "/*.toml"
		+ " -m " + getprojectfolder(METHODFOLDER) + "/method.toml"
		+ " -d 1st"
		+ " -o " + getoutputpath(RESULTFOLDER) + "/IndividualOneState")
	writeline(f,"")
	writeline(f,"#Fit all profiles using a two state exchange model (split into a several consoles to run in parrallel [remember conda environment])")
	for resi in allresidues:
			writeline(f,"chemex fit -e " + getprojectfolder(EXPERIMENTFOLDER) + "/*.toml"
				+ " -p " + getprojectfolder(PARAMETERFOLDER) + "/*.toml"
				+ " -m " + getprojectfolder(METHODFOLDER) + "/method.toml"
				+ " -o " + getoutputpath(RESULTFOLDER) + "/IndividualTwoState/" + str(resi)
				+ " --include " + str(resi))
	writeline(f,"")
	writeline(f,"#")
	writeline(f,"")
	writeline(f,"#Global two-state fitting code ([include/exclude] replace x,y,z,... with residues that should be globaly fitted):")
	writeline(f,"Go to the globalmethod file and in STEP1, select some profiles with clear minor state for initial fitting")
	writeline(f,"chemex fit -e " + getprojectfolder(EXPERIMENTFOLDER) + "/*.toml"
				+ " -p " + getprojectfolder(PARAMETERFOLDER) + "/*.toml"
				+ " -m " + getprojectfolder(METHODFOLDER) + "/globalmethod.toml"
				+ " -o " + getoutputpath(RESULTFOLDER) + "/GlobalTwoState"
				+ " --include x y z ...")
	writeline(f,"")
	writeline(f,"")
	for resi in allresidues:
			writeline(f,"chemex fit -e " + getprojectfolder(EXPERIMENTFOLDER) + "/*.toml"
				+ " -p " + getprojectfolder(PARAMETERFOLDER) + "/*.toml"
				+ " -m " + getprojectfolder(METHODFOLDER) + "/method.toml"
				+ " -d 3st"
				+ " -o " + getoutputpath(RESULTFOLDER) + "/IndividualThreeState/" + str(resi)
				+ " --include " + str(resi))
