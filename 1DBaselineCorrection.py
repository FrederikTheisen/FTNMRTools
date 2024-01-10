import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib.cm as cm
import os

SAME_WIDTH_PEAK_MODE = True
# For same width peak mode, peaks should be placed by clicking the center position
# First peak will require a second click which will define the width

args = sys.argv

if (len(args) < 2):
    print("No data file path given")
    exit()

if '-help' in args or '-h' in args:
    print("Baseline correction script, Frederik Theisen 2023")
    print("Documentation: https://github.com/FrederikTheisen/FTNMRTools")
    print()
    print()
    print("Usage: python3 1DBaselineCorrection.py <path-to-totxtexport/folder-with-totxtexports> <options>")
    print("Options:")
    print("-mode X  :   set mode, X = 0 (width mode [default]), 1 (boundary mode)")
    print("-pw X    :   set peak width, X is float")
    exit()

PATH = args[1]

BaselinePoints = {}
DataPointCount = 10
Baselines = {}
CLOSE = False
colors = {}
polydegree = 2
ppmrange = []
XAXIS = []
PEAKPOINTS = []
PEAK_WIDTH_EXPECTED = False
PEAK_WIDTH = None

if '-pw' in args:
    idx = args.index('-pw')
    PEAK_WIDTH = float(args[idx + 1])

if '-mode' in args:
    idx = args.index('-mode')
    SAME_WIDTH_PEAK_MODE = int(args[idx + 1]) == 0

def IdentifyInputData():
    if os.path.isdir(PATH): return 'dir'
    else: return 'file'

def ReadFolderData():
    global colors
    global XAXIS
    global ppmrange
    data = {}
    point = None
    minpoints = 9999999999

    files = [f for f in os.listdir(PATH) if not os.path.isdir(f) and 'Icon' not in f]

    for fn in files:
        path = PATH + "/" + fn
        spectrum = []
        with open(path) as f:
            for line in f:
                if '#' in line:
                    if '# LEFT' in line: #Get ppm range
                        dat = line.split()
                        ppmrange = [float(dat[3]),float(dat[7])]
                    if '# SIZE' in line: #Get number of data points
                        dat = line.split()
                        points = int(dat[3])
                        if points < minpoints: minpoints = points
                else: #get data point value
                    value = float(line)
                    
                    spectrum.append(value)
        data[int(fn.split('.')[0])] = spectrum

    XAXIS = list(np.linspace(ppmrange[0],ppmrange[1],minpoints))

    # Truncate data sets with too many points (number of points are +- 1 for unknown reasons)
    for d in data:
        data[d] = data[d][0:minpoints]

    # Get the 'viridis' colormap
    cmap = cm.get_cmap('viridis')

    # Generate a range of values between 0 and 1
    values = np.linspace(0, 1, len(data) + 1)

    # Map the values to colors using the colormap
    colorlist = cmap(values)

    # Remap colors to dict entries
    i = 0
    for k in data:
        colors[k] = colorlist[i]
        i += 1

    return data, minpoints

def ReadData():
    global colors
    global XAXIS
    global ppmrange
    data = {}
    points = None

    with open(PATH) as f:
        rowid = -1
        for line in f:
            if '# F2LEFT' in line:
                dat = line.split()
                ppmrange = [float(dat[3]),float(dat[7])]
            if '# NCOLS' in line:
                dat = line.split()
                points = int(dat[3])
            if "# row" in line: #Register new row
                dat = line.split()
                rowid = int(dat[-1])
                data[rowid] = []
            elif rowid > -1: #register data
                value = float(line)
                data[rowid].append(value)

    XAXIS = list(np.linspace(ppmrange[0],ppmrange[1],points))

    # Get the 'viridis' colormap
    cmap = cm.get_cmap('viridis')

    # Generate a range of values between 0 and 1
    values = np.linspace(0, 1, len(data) + 1)

    # Map the values to colors using the colormap
    colorlist = cmap(values)

    # Remap colors to dict entries
    i = 0
    for k in data:
        colors[k] = colorlist[i]
        i += 1

    return data, len(data[rowid])

def BaselineCorrect(data):
    #constants
    btn_width = 0.15
    btn_margin = 0.1

    # creating plot
    fig = plt.figure()
    ax = fig.subplots()

    plt.subplots_adjust(bottom = 0.25)

    # Add polynomium degree options
    button_ax = plt.axes([btn_margin+1*btn_width,0,btn_width,.1])
    polyplus_btn = Button(button_ax, 'Pol Degree +1', color='lightgoldenrodyellow', hovercolor='0.975')
    polyplus_btn.on_clicked(lambda event: onpolydegreebtnclick(event, data, ax, 1))
    button_ax = plt.axes([btn_margin,0,btn_width,.1])
    polyminus_btn = Button(button_ax, 'Pol Degree -1', color='lightgoldenrodyellow', hovercolor='0.975')
    polyminus_btn.on_clicked(lambda event: onpolydegreebtnclick(event, data, ax, -1))

    # Add a 'Clear' button to the plot
    button_ax = plt.axes([1-btn_margin-2*btn_width,0,btn_width,.1])
    clear_btn = Button(button_ax, 'Clear BSL Points', color="lightgoldenrodyellow", hovercolor='0.975')
    clear_btn.on_clicked(lambda event: onclearbuttonclick(event, data, ax)) # Connect the onbuttonclick function to the button
    # Add a 'Finished' button to the plot
    button_ax = plt.axes([1-btn_margin-btn_width,0,btn_width,.1])
    fin_btn = Button(button_ax, 'Subtract Baseline', color='lightgoldenrodyellow', hovercolor='0.975')
    fin_btn.on_clicked(onfinishbuttonclick) # Connect the onbuttonclick function to the button

    # Connect the onclick function to the plot
    cid = fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, data, ax))

    Draw(data, ax, mode = 0)

def PeakPicking(data):
    #constants
    btn_width = 0.15
    btn_margin = 0.1

    # creating plot
    fig = plt.figure()
    ax = fig.subplots()

    plt.subplots_adjust(bottom = 0.25)

    # Add a 'Clear' button to the plot
    button_ax = plt.axes([1-btn_margin-2*btn_width,0,btn_width,.1])
    clear_btn = Button(button_ax, 'Clear', color="lightgoldenrodyellow", hovercolor='0.975')
    clear_btn.on_clicked(lambda event: onpeakpickingclearbtnclick(event, data, ax)) # Connect the onbuttonclick function to the button
    # Add a 'Finished' button to the plot
    button_ax = plt.axes([1-btn_margin-btn_width,0,btn_width,.1])
    fin_btn = Button(button_ax, 'Finished', color='lightgoldenrodyellow', hovercolor='0.975')
    fin_btn.on_clicked(onpeakpickingfinishedclick) # Connect the onbuttonclick function to the button

    # Connect the onclick function to the plot
    cid = fig.canvas.mpl_connect('button_press_event', lambda event: onpeakpickingclick(event, data, ax))

    Draw(data, ax, mode = 1)

# Draw function
def Draw(data, ax, mode = 0):

    old_x_lim = None

    if ax.lines:
        old_x_lim = ax.get_xlim()
        old_y_lim = ax.get_ylim()

    ax.clear()
    for rowid in data:
        graph = data[rowid]
        #ax.plot(range(len(graph)), graph, color=colors[rowid])
        ax.plot(XAXIS, graph, color=colors[rowid])
        if mode == 0: #DRAW BASELINES
            if len(Baselines) > 0:
                ax.plot(XAXIS, Baselines[rowid],color=colors[rowid]) 
            
            x = []
            y = []
            
            for bp in BaselinePoints[rowid]:
                x.append(bp[0])
                y.append(bp[1])
            
            ax.scatter(x,y, color=colors[rowid])

    if mode == 1: #DRAW PEAK PICKING
        if len(PEAKPOINTS) > 0:
            print("Draw peaks")
            peakrange = ax.get_ylim()
            for peak in PEAKPOINTS:
                print(peak)
                ax.vlines(x = peak, ymin = old_y_lim[0], ymax = old_y_lim[1])

            if len(PEAKPOINTS) > 1: 
                for i in range(0,len(PEAKPOINTS) - 1,2):
                    print(i)
                    peakstart = PEAKPOINTS[i]
                    peakend = PEAKPOINTS[i+1]
                    ax.axvspan(peakstart, peakend, facecolor='b', alpha=0.2)


    if old_x_lim is not None:
        ax.set_xlim(old_x_lim)  # and restore zoom
        ax.set_ylim(old_y_lim)
    else:
        ax.invert_xaxis()

    plt.draw()
    plt.show()
    plt.pause(0.0001)

# Define a function to handle mouse clicks on the plot
def onclick(event, data, ax): 
    global Baselines
    x_pos = event.xdata
    print("Clicked at x = {}".format(x_pos))
    if event.inaxes is not None and event.dblclick:
        # Get the x position of the click
        
        print(ppmrange)
        if x_pos > ppmrange[0] or x_pos < ppmrange[1]: 
            print("ignore")
            return #not a relevant click

        AddPointsAtPosition(x_pos,data)

        Baselines = FitBaselines()

        Draw(data, ax)

# Define a functions to handle button clicks
# Baseline correction plot
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
def onpolydegreebtnclick(event, data, ax, delta):
    global polydegree
    global Baselines

    polydegree = polydegree + delta
    if polydegree < 0: polydegree = 0
    if polydegree > 20: polydegree = 20

    Baselines = FitBaselines()

    Draw(data, ax)

# Peak picking plot
def onpeakpickingclick(event, data, ax):
    x_pos = event.xdata
    print("Clicked at x = {}".format(x_pos))
    if event.inaxes is not None and event.dblclick:
        # Get the x position of the click
        
        print(ppmrange)
        if x_pos > ppmrange[0] or x_pos < ppmrange[1]: 
            print("ignore")
            return #not a relevant click

        AddPeakRangePoint(x_pos,data)

        Draw(data, ax, mode = 1)
def onpeakpickingclearbtnclick(event, data, ax):
    PEAKPOINTS.clear()
    print(len(PEAKPOINTS))

    Draw(data, ax, mode = 1)
def onpeakpickingfinishedclick(event):
    print("done")
    plt.close()

def AddPointsAtPosition(position,data):
    print(GetAxisIndexFromPosition(position))
    for rowid in data:
        BaselinePoints[rowid].append([position,data[rowid][GetAxisIndexFromPosition(position)]])

def GetAxisIndexFromPosition(position):
    for i in range(len(XAXIS)):
        pos = XAXIS[i]
        if pos < position: return i

def FitBaselines():
    baselines = {}

    for rowid in BaselinePoints:
        if len(BaselinePoints[rowid]) < polydegree + 1: return baselines

        x = []
        y = []
        for point in BaselinePoints[rowid]:
            x.append(point[0])
            y.append(point[1])

        fit = FitBaseline(x,y)

        baselines[rowid] = []
        for i in XAXIS:
            baselines[rowid].append(np.polyval(fit,i))

    return baselines

def FitBaseline(x,y):
    return np.polyfit(x,y,polydegree)

def SubtractBaseline(data):
    baselinecorrected = {}

    for rowid in data:
        bsl = Baselines[rowid]
        dat = data[rowid]
        baselinecorrected[rowid] = []
        for i in range(len(bsl)):
            baselinecorrected[rowid].append(dat[i]-bsl[i])

    return baselinecorrected

def AddPeakRangePoint(position, data):
    global PEAK_WIDTH
    global PEAK_WIDTH_EXPECTED
    axis_idx = GetAxisIndexFromPosition(position)

    if SAME_WIDTH_PEAK_MODE:
        print("peak width mode")
        if PEAK_WIDTH is None and not PEAK_WIDTH_EXPECTED: # No peak width defined, place peak marker
            print("set peak center")
            PEAK_WIDTH_EXPECTED = True
            PEAKPOINTS.append(position)
            return
        elif PEAK_WIDTH is None and PEAK_WIDTH_EXPECTED: # Peak width setup, clear existing PEAKPOINTS
            print("set peak width")
            PEAK_WIDTH_EXPECTED = False
            PEAK_WIDTH = abs(PEAKPOINTS[0] - position) # Calc peak width
            print("PEAK_WIDTH = " + str(PEAK_WIDTH))
            print("-pw " + str(PEAK_WIDTH))
            position = PEAKPOINTS[0] # tmp save first click
            PEAKPOINTS.clear() # clear points

        peak_start = position - PEAK_WIDTH # setup peak
        peak_end = position + PEAK_WIDTH

        PEAKPOINTS.append(peak_end)
        PEAKPOINTS.append(peak_start)
            
    else:
        print("peak border mode")
        PEAKPOINTS.append(position)

### EXPORT FUNCTIONS ###

def PrintData(dat):
    header = "ppm "
    for rowid in dat:
        header += str(rowid) + " "
    with open("baselined.txt","w+") as f:
        f.write(header + "\n")
        for i in range(DataPointCount):
            out = str(XAXIS[i]) + " "
            for rowid in dat:
                row = dat[rowid]
                out += str(row[i]) + " "
            f.write(out.strip() + "\n")

def ExportPeakVolumes(data):
    volumes = {}
    peakcount = 0

    for rowid in data: volumes[rowid] = []

    for i in range(0,len(PEAKPOINTS) - 1,2):
        peakstart = PEAKPOINTS[i]
        peakend = PEAKPOINTS[i+1]
        idx_start = GetAxisIndexFromPosition(peakstart)
        idx_end = GetAxisIndexFromPosition(peakend)
        peakcount += 1

        for rowid in data:
            volume = sum(data[rowid][idx_start:idx_end])
            volumes[rowid].append(volume)

    header = "peak "
    for rowid in data:
        header += str(rowid) + " "

    with open("peaks.txt","w+") as f:
        f.write(header + "\n")
        for i in range(peakcount):
            out = str(i) + " "
            for rowid in volumes:
                vs = volumes[rowid]
                out += str(vs[i]) + " "
            f.write(out.strip() + "\n")

### MAIN LOOP ###
    
def Main():
    global DataPointCount

    if IdentifyInputData() == 'dir': data,DataPointCount = ReadFolderData()
    else: data,DataPointCount = ReadData()

    for rowid in data: BaselinePoints[rowid] = []

    BaselineCorrect(data)

    print("Subtracting baselines...")
    corr = SubtractBaseline(data)
    plt.close()
    print("Saving baseline corrected data...")
    PrintData(corr)
    print("Done")

    PeakPicking(corr)

    ExportPeakVolumes(corr)

Main()
