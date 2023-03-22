import sys
from os import listdir
from os.path import isfile,join
import tomli

PATH = sys.argv[1]

readdata = input("Read Data Files [y/n]: ").lower() == 'y'
resi = input("Residue index: ")
ext = '.fit' 
if readdata: ext = '.exp'

xvals = []

data = {}

files = [f for f in listdir(PATH) if isfile(join(PATH, f)) and f.endswith(ext)]

print(files)


for file in files:
	with open(PATH + "/" + file) as f:
		lines = f.readlines()
		datn = ""
		values = {}
		for line in lines:
			if len(line) < 2: continue
			if line[0] == '[':
				if datn is not "": 
					data[datn] = values
					values = {}
				datn = line[1:-2]
			elif line[0] != '#':
				dat = line.split()
				x = float(dat[0])
				if x not in xvals: xvals.append(x)
				values[x] = float(dat[1])
		data[datn] = values

	header = "ppm "
	for datn in data:
		header += datn + " "
	print(header)

	for x in sorted(xvals):
		line = [str(x)]
		for datn in data:
			dat = data[datn]
			if x in dat:
				line.append(str(dat[x]))
			else: line.append(" ")
		print(' '.join(line))



