import os
import time

# Function to get the data of the experiment recording
def GetTime(path):
    # Set aqcus file path
    acqus_path = path + "/acqus"
    # Read the date and time from the acqu file
    with open(acqus_path, 'r') as f:
        # Loop through all lines in the file until a line stars with "##$DATA=", then return the value after the equals sign
        for line in f:
            if line.startswith('##$DATE='):
                return int(line.split('=')[1].strip())

# Define the path to the directory to search
dir_path = os.getcwd()

# Define a dictionary to store the 'fid' creation time for each folder
fid_times = {}

# Loop through all directories in the directory
for folder_name in os.listdir(dir_path):
    # Check if folder/file name is a number and if the number is larger than 1000 (time traces use expno > 1000
    if folder_name.isdigit() and int(folder_name) > 1000:
        # Construct path to folder
        folder_path = os.path.join(dir_path, folder_name)
        # Get exp time from folder
        fid_times[int(folder_name)] = GetTime(folder_path)

# Print how many FID times were found (debug info)
print("Logged FID times: " + str(len(fid_times)))

# Set first start experiment (series start with expno x001)
refn = 1001

# while the selected series (refn) exists, do this:
while refn in fid_times:
    # Get time for first exp of series  
    reft = fid_times[refn]
    print("Exp #" + str(refn) + ", time: " + str(reft))
    # Loop through all exp in the current series (x001 to x999)
    for i in range(refn,refn+999):
        if i in fid_times: # check if FID time logged 
            # get the current exp time difference to the first exp
            offset_time = fid_times[i] - reft
            print(i,offset_time)
        else: break #if expno not found, then stop this loop (THIS LINE MAKES THE SCRIPT INCOMPATIBLE WITH MISSING EXPNOs)
    # Go to next series
    refn += 1000
