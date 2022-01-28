from math import atan2, degrees, tan, pi, sqrt, log10, radians
import numpy as np

# GOMspace NanoCam C1U with MT9T031 CMOS sensor specs

h = 6.55 #mm
v = 4.92 #mm
size = h*v
h_pixels = 2048
v_pixels = 1536
total_pixels = h_pixels*v_pixels
pixel_size = 3.2*3.2 #um^2
focal_length = 8 #mm f/1.4
frame_rate = 12 #fps 
output = 10 #bit
output_gain = 25 #e-/LSB
read_noise = 6 #e-RMS
dark_current = 100 #e-/pixel/s
quantum_efficiency = .25

# FOV, swath, resolution calculations

altitude_min = 500 #km
altitude_max = 3200 #km

FOV_h = 2*atan2(.5*h,focal_length) #radians
FOV_v = 2*atan2(.5*v,focal_length) #radians

swath_h_min = 2*altitude_min*tan(FOV_h/2) #km
swath_v_min = 2*altitude_min*tan(FOV_v/2) #km

swath_h_max = 2*altitude_max*tan(FOV_h/2) #km
swath_v_max = 2*altitude_max*tan(FOV_v/2) #km

resolution_h_min = swath_h_min*1e3/h_pixels #m
resolution_v_min = swath_v_min*1e3/v_pixels #m

resolution_h_max = swath_h_max*1e3/h_pixels #m
resolution_v_max = swath_v_max*1e3/v_pixels #m

# Image compression

def compression(alt, h_bin_size, v_bin_size, time, h_pixels=2048, v_pixels=1536, focal_length=8, h=6.55, v=4.92, output=10, frame_rate=12):
    
    FOV_h = 2*np.arctan2(.5*h,focal_length)
    FOV_v = 2*np.arctan2(.5*v,focal_length)
    
    swath_h = 2*alt*np.tan(FOV_h/2)
    swath_v = 2*alt*np.tan(FOV_v/2)
    swaths = np.array([swath_h,swath_v]).T
    
    res_conversion = np.array([1e3/h_pixels,1e3/v_pixels])
    res = swaths*res_conversion
    
    binned_res = np.zeros(np.shape(res))
    h_bins = h_pixels/h_bin_size
    binned_res[:,0] = res[:,0]*h_bin_size
    v_bins = v_pixels/v_bin_size
    binned_res[:,1] = res[:,1]*v_bin_size
    
    data = output*frame_rate*h_pixels*v_pixels*time*60/8e9
    data_rate = output*frame_rate*h_pixels*v_pixels
    data_binned = output*frame_rate*h_bins*v_bins*time*60/8e9
    
    data_compressed = output*(1/20)*h_bins*v_bins*time*60/8e9
    
    return (res, binned_res, data, data_binned, data_rate, data_compressed)

altitudes = np.arange(500,3300,100)
data_compression = compression(altitudes, 8, 4, 360)

# Sensitivity Calculations

## SIGNAL
optical_E = 3.6e6 #J (flash)
h = 6.62607004e-34 #m^2 kg/s
l = 777.4 %nm
E_photon = h*(3e8)/(l*1e-9) #J
photons_per_flash = optical_E/E_photon #photons

flux_at_SC_min = photons_per_flash/(4*pi*(altitude_min**2)) #photons/km^2
flux_at_SC_max = photons_per_flash/(4*pi*(altitude_max**2)) #photons/km^2

photons_at_detector_min = flux_at_SC_min*(size*1e-12) #photons (convert size from mm^2 to km^2)
photons_at_detector_max = flux_at_SC_max*(size*1e-12) #photons (convert size from mm^2 to km^2)

signal_min = photons_at_detector_min*quantum_efficiency
signal_max = photons_at_detector_max*quantum_efficiency

## BACKGROUND
background_at_SC_min = 9.88e14 #at min altitude
background_at_SC_max = 2.4e13 #at max altitude

background_at_detector_min = background_at_SC_min*(size*1e-12)
background_at_detector_max = background_at_SC_max*(size*1e-12)

background_noise_min = background_at_detector_min*quantum_efficiency
background_noise_max = background_at_detector_max*quantum_efficiency

## SNR
total_noise_min = sqrt((output_gain*2000)**2 + read_noise**2 + dark_current**2 + background_noise_min**2)
total_noise_max = sqrt((output_gain*2000)**2 + read_noise**2 + dark_current**2 + background_noise_max**2)

SNR_min = signal_min/total_noise_min
SNR_max = signal_max/total_noise_max

SNR_dB_min = 20*log10(SNR_min)
SNR_dB_max = 20*log10(SNR_max)
