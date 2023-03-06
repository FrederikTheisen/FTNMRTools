import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib.cm as cm

args = sys.argv

if (len(args) < 2):
	print("No data file path given")
	exit()

PATH = args[1]

BaselinePoints = {}
DataPointCount = 10
Baselines = {}
CLOSE = False
colors = []

def ReadData():
	global colors
	data = {}

	with open(PATH) as f:
		rowid = -1
		for line in f:
			if "# row" in line: #Register new row
				dat = line.split()
				rowid = int(dat[-1])
				data[rowid] = []
			elif rowid > -1: #register data
				value = float(line)
				data[rowid].append(value)

	# Get the 'viridis' colormap
	cmap = cm.get_cmap('viridis')

	# Generate a range of values between 0 and 1
	values = np.linspace(0, 1, len(data) + 1)

	# Map the values to colors using the colormap
	colors = cmap(values)

	return data, len(data[rowid])

def PrintData(dat):
	header = "x "
	for rowid in dat:
		header += str(rowid) + " "
	print(header)

	for i in range(DataPointCount):
		out = str(i) + " "
		for rowid in dat:
			row = dat[rowid]
			out += str(row[i]) + " "
		print(out.strip())

def BaselineCorrect(data):
	
	#constants
	btn_width = 0.15
	btn_margin = 0.1

	# creating plot
	fig = plt.figure()
	ax = fig.subplots()
	plt.subplots_adjust(bottom = 0.25)

	# Add polynomium degree options
	# button_ax = plt.axes([1-btn_margin-2*btn_width,0,btn_width,.1])
	# clear_btn = Button(button_ax, 'Pol Degree +1', color='lightgoldenrodyellow', hovercolor='0.975')
	# clear_btn.on_clicked(lambda event: onclearbuttonclick(event, data, ax))
	# button_ax = plt.axes([1-btn_margin-2*btn_width,0,btn_width,.1])
	# clear_btn = Button(button_ax, 'Pol Degree -1', color='lightgoldenrodyellow', hovercolor='0.975')
	# clear_btn.on_clicked(lambda event: onclearbuttonclick(event, data, ax))

	# Add a 'Clear' button to the plot
	button_ax = plt.axes([1-btn_margin-2*btn_width,0,btn_width,.1])
	clear_btn = Button(button_ax, 'Clear', color="lightgoldenrodyellow", hovercolor='0.975')
	clear_btn.on_clicked(lambda event: onclearbuttonclick(event, data, ax)) # Connect the onbuttonclick function to the button
	# Add a 'Finished' button to the plot
	button_ax = plt.axes([1-btn_margin-btn_width,0,btn_width,.1])
	fin_btn = Button(button_ax, 'Finished', color='lightgoldenrodyellow', hovercolor='0.975')
	fin_btn.on_clicked(onfinishbuttonclick) # Connect the onbuttonclick function to the button

	# Connect the onclick function to the plot
	cid = fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, data, ax))

	Draw(data, ax)


	# Show the plot and wait for user interaction
	#plt.show()

def Draw(data, ax):

	old_x_lim = None

	if ax.lines:
		old_x_lim = ax.get_xlim()
		old_y_lim = ax.get_ylim()

	ax.clear()
	for rowid in data:
		graph = data[rowid]
		ax.plot(range(len(graph)), graph, color=colors[rowid])
		if len(Baselines) > 0:
			print("Draw baselines")
			ax.plot(range(len(graph)), Baselines[rowid],color=colors[rowid])
		x = []
		y = []
		for bp in BaselinePoints[rowid]:
			x.append(bp[0])
			y.append(bp[1])
		ax.scatter(x,y, color=colors[rowid])


	if old_x_lim is not None:
		ax.set_xlim(old_x_lim)  # and restore zoom
		ax.set_ylim(old_y_lim)

	plt.draw()
	plt.show()
	plt.pause(0.0001)
	plt.clf()

def onclick(event, data, ax): # Define a function to handle mouse clicks on the plot
	global Baselines

	if event.inaxes is not None and event.dblclick:
	    # Get the x position of the click
	    x_pos = event.xdata
	    if x_pos < 1: return #not a relevant click
	    print("Clicked at x = {}".format(x_pos))

	    AddPointsAtPosition(x_pos,data)

	    Baselines = FitBaselines()

	    Draw(data, ax)

# Define a function to handle button clicks
def onfinishbuttonclick(event):
    #global CLOSE
    print("Finished button clicked")
    #CLOSE = True
    plt.close()
def onclearbuttonclick(event, data, ax):
    global BaselinePoints
    global Baselines
    print("Clear button clicked")
    for rowid in data: BaselinePoints[rowid] = []
    Baselines = {}

    Draw(data, ax)

def AddPointsAtPosition(position,data):
	position = int(position)

	for rowid in data:
		BaselinePoints[rowid].append([position,data[rowid][position]])

def FitBaselines():
	baselines = {}

	for rowid in BaselinePoints:
		if len(BaselinePoints[rowid]) < 3: return baselines

		x = []
		y = []
		for point in BaselinePoints[rowid]:
			x.append(point[0])
			y.append(point[1])

		fit = FitBaseline(x,y)

		baselines[rowid] = []
		for i in range(DataPointCount):
			baselines[rowid].append(np.polyval(fit,i))

	return baselines

def FitBaseline(x,y):
	return np.polyfit(x,y,2)

def SubtractBaseline(data):
	baselinecorrected = {}

	for rowid in data:
		bsl = Baselines[rowid]
		dat = data[rowid]
		baselinecorrected[rowid] = []
		for i in range(len(bsl)):
			baselinecorrected[rowid].append(dat[i]-bsl[i])

	return baselinecorrected


def Main():
	global DataPointCount

	data,DataPointCount = ReadData()

	for rowid in data: BaselinePoints[rowid] = []

	BaselineCorrect(data)

	corr = SubtractBaseline(data)

	PrintData(corr)

Main()
