import sys
import os
from os import listdir
from os.path import isfile, join
import math
import tomli
import matplotlib.pyplot as plt

PATH = sys.argv[1]

def getpathforfolder(folder):
	return PATH + "/" + folder

def fixtomltrash(tomlfiledata):
	newfile = []
	for line in tomlfiledata:
		if '±' in line:
			data = line.split()
			line = data[0] + " = " + "[" + data[2] + ", " + data[4][1:] + "]"
		newfile.append(line)
	return '\n'.join(newfile)

def readtomlkey(data, key):
	var = list(data[key].values())[0]
	if type(var) is float: var = [var, 0]
	return var

def readtomllistvalue(value):
	if type(value) is float: 
		return [value, 0]
	else: 
		return value

resfolders = [f for f in listdir(PATH) if not isfile(join(PATH,f))]

minmax = [999999,-999999]
results = {}

for folder in sorted(resfolders):
	if int(folder) > minmax[1]: minmax[1] = int(folder)
	if int(folder) < minmax[0]: minmax[0] = int(folder)
	path = getpathforfolder(folder)

	maxstep = 0
	for d in listdir(path): #Handle multistep processing
		if 'STEP' in d: 
			if int(d[4:]) > maxstep: maxstep = int(d[4:])
	if maxstep > 0: path += "/STEP" + str(maxstep)

	params = path + '/Parameters/fitted.toml'
	plot = join(path,'/Plot')
	stats = path + "/statistics.toml"
	result = {}
	if True:
		with open(stats) as f:
			lines = f.readlines()
			result["rchi"] = float(lines[3].split()[2])
		with open(params) as f:
			data = tomli.loads(fixtomltrash(f.readlines()))
			globaldata = data['GLOBAL']
			for key in globaldata.keys(): result[key] = readtomllistvalue(globaldata[key])

			result['CS_A'] = readtomlkey(data, "CS_A")
			result['DW_AB'] = readtomlkey(data, "DW_AB")
			result['CS_B'] = [result['CS_A'][0] + result['DW_AB'][0], math.sqrt(result['CS_A'][1]**2 + result['DW_AB'][1]**2)]
			if 'DW_AC' in data:
				result['DW_AC'] = readtomlkey(data, "DW_AC")
				result['CS_C'] = [result['CS_A'][0] + result['DW_AC'][0], math.sqrt(result['CS_A'][1]**2 + result['DW_AC'][1]**2)]
		results[int(folder)] = result


for someresidueindex in results.keys():
	columns = ['peak']
	for key in results[someresidueindex].keys():
		val = results[someresidueindex][key]
		columns.append(key)
		if type(val) is not float: columns.append(key + "_ERROR")
	print(' '.join(columns))
	break

for resi in range(minmax[0],minmax[1]+1):
	row = [str(resi)]
	if resi in results.keys():
		for key in results[resi]:
			value = results[resi][key]
			if type(value) is float: 
				row.append(str(results[resi][key]))
			else: 
				row.append(str(results[resi][key][0]))
				row.append(str(results[resi][key][1]))

	print(' '.join(row))


exit()
### CLUSTERING ###
fig = plt.figure()
ax = fig.add_subplot(2, 1, 1)

for resi in results:
	#if resi not in [263,265,267,268]: continue
	r = results[resi]
	x = [r[k][0] for k in r if 'KEX' in k and 'ERROR' not in k]
	y = [r[k][0] for k in r if 'P' in k and len(k) == 2]
	ax.plot(x, y, '-', label = "resi " + str(resi))
	x = [r['KEX_AB'][0]]
	y = [r['PB'][0]]
	ax.plot(x, y, 'o', label = "resi " + str(resi))

ax.set_xscale('log')
#ax.set_yscale('log')
plt.legend()
plt.show()
	
