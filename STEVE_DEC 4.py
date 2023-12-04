
# THIS CODE WILL PRODUCE ION TEMPERATURE FIGURES FOR LONG-PULSE DATAFILES 
# DOWNLOADED FROM "data.amisr.com" AND THEN STRIPPED USING "Code-for-Andy-1.py"


# IMPORT THE REQUIRED LIBRARIES (DO NOT TOUCH, UNLESS ALTERING THE CODE)
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import math
import datetime
import time
import sys
from scipy.signal import find_peaks
import os


# CUSTOMIZE THE PATH TO THE DATAFILE RELATIVE TO THE CODE
#path = "Andy PFISR/"

# INSERT THE NAME OF THE h5 DATAFILE YOU HAVE DOWNLOADED AND WOULD LIKE TO PLOT
#filename = "TIspike20180606.003_lp_5min-fitcal.h5"

# THIS OPENS THE DATA FILE OF INTEREST
#f  = h5.File(path+filename, 'r')


# Absolute path to the directory containing the data file
path = r"C:\Users\armee\OneDrive\Documents\Research work\Andy PFISR"

# Name of the h5 data file you want to plot
filename = "TIspike20180606.003_lp_5min-fitcal.h5"

# Construct the full file path
file_path = os.path.join(path, filename)

# Open the data file
f = h5.File(file_path, 'r')


# BEAM PARAMETERS (E.G. ELEVATION ANGLE, AZIMUTH, ETC)
beamcodes = f['BeamCodes']
beamcodesdata = beamcodes[:,:]
beamcodesdata = beamcodesdata.astype(float)    

# EPOCH TIME
epoch = f['UnixTime']
epochData1 = epoch[:,:]
epochData1 = epochData1.astype(float)    
epochData=[]
for i in range(len(epochData1)): 
   epochData.extend([np.mean([epochData1[i,0],epochData1[i,1]])]) 

# ALTITUDE
alt = f['Altitude']
altdata = alt[:,:]/1000.#CONVERTING FROM METERS TO KILOMETERS
altdata = altdata.astype(float)    

# PLASMA DENSITY [m-3]
ne = f['Ne']
nedata = ne[:,:,:]
nedata = nedata.astype(float)    

# ERROR IN PLASMA DENSITY [m-3]
dne = f['dNe']
dnedata = dne[:,:,:]
dnedata = dnedata.astype(float)    

# THE ION SPECIES IS SET
# 0 = O+, THE DOMINANT F-REGION ION. SOME DATAFILES OFFER, 1 = O2+, 2 = NO+ , 
# 3 = N2+ AND , 4 = N+, BUT I STRONGLY SUGGEST YOU DO NOT TOUCH THIS
ion = 0

# ION TEMPERATURE [K]
ti = f['Fits']
tidata = ti[:,:,:,ion,1]
tidata = tidata.astype(float)    

# ERROR IN ION TEMPERATURE [K]
dti = f['Errors']
dtidata = dti[:,:,:,ion,1]
dtidata = dtidata.astype(float)    

# ELECTRON TEMPERATURE [K]
te = f['Fits']
tedata = te[:,:,:,-1,1]  
tedata = tedata.astype(float)    

# ERROR IN ELECTRON TEMPERATURE [K]
dte = f['Errors']
dtedata = dte[:,:,:,-1,1]
dtedata = dtedata.astype(float)    

# LINE-OF-SIGHT ION VELOCITY [m/s]
losv = f['Fits']
losvdata = losv[:,:,:,ion,3]
losvdata = losvdata.astype(float)    

# CORRECTED GEOMAGNETIC LATITUDE [deg]
lat = f['MagneticLatitude']
latdata = lat[:,:]
latdata = latdata.astype(float)       

# CORRECTED GEOMAGNETIC LONGITUDE [deg]
long = f['MagneticLongitude']
longdata = long[:,:]
longdata = longdata.astype(float)  

# THIS CLOSES THE DATA FILE
f.close() 

# THIS CONVERTS EPOCH TIME TO UNIVERSAL TIME IN DATETIME
ut = [datetime.datetime(int(time.gmtime((epochData1[t,0]+epochData1[t,1])/2)[0]), # YEAR
                       int(time.gmtime((epochData1[t,0]+epochData1[t,1])/2)[1]), # MONTH
                       int(time.gmtime((epochData1[t,0]+epochData1[t,1])/2)[2]), # DAY
                       int(time.gmtime((epochData1[t,0]+epochData1[t,1])/2)[3]), # HOUR
                       int(time.gmtime((epochData1[t,0]+epochData1[t,1])/2)[4]), # MINUTE
                       int(time.gmtime((epochData1[t,0]+epochData1[t,1])/2)[5])) # SECOND
                       for t in range(len(epochData1))]




# =============================================================================
# # Define a function to identify peaks and troughs
# def find_multiple_peaks(Ti):
#         # Smooth the ion temperature data for better peak detection
#         window_size = 5
#         smoothed_Ti = np.convolve(Ti, np.ones(window_size)/window_size, mode='same')
# 
#         # Calculate threshold for peak detection (adjust as needed)
#         threshold = np.nanmean(smoothed_Ti) + 3 * np.nanstd(smoothed_Ti)
# 
#         # Find peaks in the smoothed data using SciPy's find_peaks function
#         peaks, _ = find_peaks(smoothed_Ti, height=threshold, distance=20)
# =============================================================================

# =============================================================================
#         return peaks
# =============================================================================
# FIXME: identify multiple peaks
def find_multiple_peaks(Ti, min_distance=20, relative_height=0.02):
    # Normalize the data to have values between 0 and 1
    Ti_normalized = (Ti - np.nanmin(Ti)) / (np.nanmax(Ti) - np.nanmin(Ti))

    # Find peaks in the normalized data
    peaks, _ = find_peaks(Ti_normalized, distance=min_distance, height=relative_height)

    return peaks

# Initialize lists to store ion temperatures and peak indices for each altitude
Ti_210, Ti_240, Ti_270 = [], [], []
peaks_210, peaks_240, peaks_270 = [], [], []

# Loop through the beams
for i in range(len(altdata)):
    if beamcodesdata[i, 0] == 64157.0:
        # Process for each target altitude
        for target_altitude, Ti_all, peaks_all in zip([210, 240, 270], 
                                                      [Ti_210, Ti_240, Ti_270], 
                                                      [peaks_210, peaks_240, peaks_270]):
            
            # Find ion temperature closest to the target altitude
            altitude_diff = np.abs(altdata[i, :] - target_altitude)
            closest_altitude_index = np.nanargmin(altitude_diff)

            # Extract ion temperature at the closest altitude
            Ti = tidata[:, i, closest_altitude_index]
            Ti_all.append(Ti)

            # Find peaks
            peaks = find_multiple_peaks(Ti)
            peaks_all.append(peaks)

# Aggregate data for plotting
# ... [rest of your aggregation and plotting code]


# Plot ion temperature data with highlighted peaks
# Plot ion temperature data with highlighted peaks
plt.figure(figsize=(14, 7))  # Adjust figure size as needed

Ti_210 = np.squeeze(Ti_210)
Ti_240 = np.squeeze(Ti_240)
Ti_270 = np.squeeze(Ti_270)


# Calculate the min and max temperatures across all datasets


# Plot ion temperatures for each altitude
plt.plot(ut, Ti_210, color="red", label='Ion Temp - 210 km')
plt.plot(ut, Ti_240, color="green", label='Ion Temp - 240 km')
plt.plot(ut, Ti_270, color="blue", label='Ion Temp - 270 km')

# Highlight the peaks
plt.scatter(np.array(ut)[peaks_210], np.array(Ti_210)[peaks_210], color="red", marker="o", label="Peaks - 210 km")
plt.scatter(np.array(ut)[peaks_240], np.array(Ti_240)[peaks_240], color="green", marker="o", label="Peaks - 240 km")
plt.scatter(np.array(ut)[peaks_270], np.array(Ti_270)[peaks_270], color="blue", marker="o", label="Peaks - 270 km")

# Set plot title and labels
plt.title('Ion Temperature Peaks Identified - ' + filename)
plt.xlabel('Universal Time')
plt.ylabel('Ion Temperature [K]')

# Set y-axis limits to focus on the range of interest
# Adjust these values based on your data range or desired focus
plt.ylim([0, 50000])  # Replace with  min and max values

# Add grid and legend
plt.grid(True)
plt.legend()

# Display the plot
plt.show()

# Optionally, save the plot to a file
# plt.savefig('your_plot_filename.png')  # Replace with your desired file path and name

# Clear the current figure after displaying/saving
plt.cla()
plt.clf()
plt.close()

# # Loop through the beams
# for i in range(len(altdata)):
#     if beamcodesdata[i, 0] == 64157.0:
#         # Find ion temperature closest to 210 km
#         target_altitude = 210
#         altitude_diff = np.abs(altdata[i, :] - target_altitude)
#         closest_altitude_index = np.nanargmin(altitude_diff)
#         # Extract ion temperature at the closest altitude
#         Ti_210 = tidata[:, i, closest_altitude_index]
#         peaks_210 = find_multiple_peaks(Ti_210)
# 
#         # Find ion temperature closest to 240 km
#         target_altitude = 240
#         altitude_diff = np.abs(altdata[i, :] - target_altitude)
#         closest_altitude_index = np.nanargmin(altitude_diff)
#         # Extract ion temperature at the closest altitude
#         Ti_240 = tidata[:, i, closest_altitude_index]
#         peaks_240 = find_multiple_peaks(Ti_240)
# 
#         # Find ion temperature closest to 270 km
#         target_altitude = 270
#         altitude_diff = np.abs(altdata[i, :] - target_altitude)
#         closest_altitude_index = np.nanargmin(altitude_diff)
#         # Extract ion temperature at the closest altitude
#         Ti_270 = tidata[:, i, closest_altitude_index]
#         peaks_270 = find_multiple_peaks(Ti_270)
# # FIXME: here, put in code that will compare the peaks for the three altitudes and remove peaks that are not consistent
# 
# # Loop through the beams
# for i in range(len(altdata)):
#     if beamcodesdata[i, 0] == 64157.0:
#         # Process for each target altitude
#         for target_altitude in [210, 240, 270]:
#             # Find ion temperature closest to the target altitude
#             altitude_diff = np.abs(altdata[i, :] - target_altitude)
#             closest_altitude_index = np.nanargmin(altitude_diff)
# 
#             # Extract ion temperature at the closest altitude
#             Ti = tidata[:, i, closest_altitude_index]
#             peaks = find_multiple_peaks(Ti)
# 
#             # Further processing or plotting can be done here with 'peaks'
#             # For example, plotting or analyzing these peaks.
# 
# 
# # Construct the plot file name
# #if len(sys.argv) > 1:
#     #plotfile = 'Andy PFISR/230{} - {}.png'.format(filename.split("_")[0], sys.argv[1])
# #else:
#     #plotfile = 'Andy PFISR/230{} - NoArg.png'.format(filename.split("_")[0])
# 
# # ... [rest of your code]
# 
# # Plot ion temperature data with highlighted peaks
# plt.figure(figsize=(14, 7))
# plt.plot(ut, Ti_210, color="red", label='Raw Ion Temp - 210 km')
# plt.plot(ut, Ti_240, color="green", label='Raw Ion Temp - 240 km')
# plt.plot(ut, Ti_270, color="blue", label='Raw Ion Temp - 270 km')
# plt.scatter(np.array(ut)[peaks_210], np.array(Ti_210)[peaks_210], color="red", marker="o", label="Peaks - 210 km")
# plt.scatter(np.array(ut)[peaks_240], np.array(Ti_240)[peaks_240], color="green", marker="o", label="Peaks - 240 km")
# plt.scatter(np.array(ut)[peaks_270], np.array(Ti_270)[peaks_270], color="blue", marker="o", label="Peaks - 270 km")
# plt.title('Ti Spike Identified_' + filename)
# plt.grid(True)
# plt.xlabel('Universal Time')
# plt.ylabel('Ion Temperature [K]')
# plt.legend()
# 
# # Display or save the plot
# plt.show()  # Use this to display the plot
# # plt.savefig(plotfile)  # Uncomment to save the plot
# 
# # Clear the current figure
# plt.cla()
# plt.clf()
# plt.close()
# 
# =============================================================================
