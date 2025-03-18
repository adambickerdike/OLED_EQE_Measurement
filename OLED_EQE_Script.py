# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 12:52:54 2025

@author: bicke
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OLED Characterization System using SCPI
Enhanced version with improved efficiency and organization

Original author: Adam Bickerdike
Refactored: 2025-03-18

This script implements a comprehensive measurement system for organic light-emitting diodes (OLEDs).
It interfaces with a Keysight SMU and Ocean Optics spectrometer to characterize:
- Current-voltage characteristics
- Luminance and spectral output
- Quantum efficiency
- Color coordinates
- Device mobility and threshold voltage
"""

import sys
import time
import signal
import traceback
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import constants
from scipy.interpolate import interp1d
from astropy.convolution import convolve, Box1DKernel
import pyvisa as visa
import pwlf
import seabreeze
from seabreeze.spectrometers import Spectrometer
from contextlib import contextmanager

# Configure matplotlib for better plot appearance
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

# Initialize seabreeze
seabreeze.use('pyseabreeze')

class OLEDMeasurement:
    def __init__(self):
        # Device parameters
        self.dev_width = 0.025  # cm
        self.dev_length = 0.01  # cm
        self.device_area = self.dev_width * self.dev_length  # cm²
        self.device_area_m = self.device_area * 1e-4  # m²
        
        # Measurement settings
        self.start_voltage = 0.0  # V
        self.stop_voltage = 10.0  # V
        self.step_size = 0.5  # V
        self.calibration_int_time_ms = 5.0  # ms
        self.spec_integration_time = 1000.0  # ms
        self.current_measure_time = 0.5  # seconds
        self.compliance_current = 0.1  # A
        
        # Output file paths
        self.filepath = 'D:\\Integrating_Sphere_Results\\06_09_2024\\'
        self.filename_results = 'BSBCz_uOLED_250um_1_1.csv'
        self.filename_wavelength = 'BSBCz_uOLED_250um_Spectra_1_1.csv'
        self.filename_spectra = 'BSBCz_uOLED_250um_Spectra_1.csv'
        
        # Spectrometer settings
        self.boxcar_width = 5
        self.spectra_repeats = 3
        self.optimal_rms = 500
        
        # Initialize data arrays
        self.init_data_arrays()
        
        # Setup interrupt handling
        self.interrupted = False
        signal.signal(signal.SIGINT, self.signal_handler)
        
    def init_data_arrays(self):
        """Initialize empty arrays for data collection"""
        self.voltage = []
        self.current = []
        self.current_density = []
        self.brightness = []
        self.current_efficiency = []
        self.power_efficiency = []
        self.photocurrent = []
        self.eqe = []
        self.cie_x = []
        self.cie_y = []
        self.spectra = []
        self.interp_spectra = []
        self.corrected_spectra = []
        self.power_spectra = []
        self.lum_list_nm = []
        self.luminance = []
        
        # Calibration data
        self.cie_1nm_x = []
        self.cie_1nm_y = []
        self.cal_x = []
        self.cal_y = []
        self.cie_xyz_l = []
        self.cie_xyz_x = []
        self.cie_xyz_y = []
        self.cie_xyz_z = []
        
    def signal_handler(self, signal, frame):
        """Handle Ctrl+C interruption"""
        self.interrupted = True
        print("\nMeasurement interrupted by user. Finishing...")
    
    def load_calibration_data(self):
        """Load calibration data from files"""
        # CIE data
        cie_data = pd.read_csv('CIE_1nm_NIR_test.csv')
        self.cie_1nm_x = cie_data['CIE_x'].values[20:421]
        self.cie_1nm_y = cie_data['CIE_y'].values[20:421]
        
        # Calibration curve
        cal_data = pd.read_csv('1000um_Calibration_Curve_27_03_2022.csv')
        self.cal_x = cal_data['Wavelength'].values
        self.cal_y = cal_data['Value'].values
        
        # CIE xyz data
        cie_xyz = pd.read_csv('CIE_xyz.csv')
        self.cie_xyz_l = cie_xyz['l'].values
        self.cie_xyz_x = cie_xyz['x'].values
        self.cie_xyz_y = cie_xyz['y'].values
        self.cie_xyz_z = cie_xyz['z'].values
    
    def connect_instruments(self):
        """Connect to instruments via VISA"""
        self.rm = visa.ResourceManager()
        
        # Connect to the spectrometer
        print('Connecting to spectrometer...')
        try:
            self.spec = Spectrometer.from_first_available()
            print(f'Spectrometer connected: {self.spec}')
        except Exception as e:
            print(f"Error connecting to spectrometer: {e}")
            sys.exit(1)
            
        # Connect to Keysight SMU
        print('Connecting to SMUs...')
        try:
            self.smu = self.rm.open_resource('GPIB0::24::INSTR')
            self.smu.timeout = 250000
            print(f'SMU connected: {self.smu.query("*IDN?")}')
        except Exception as e:
            print(f"Error connecting to SMU: {e}")
            sys.exit(1)
    
    def initialize_smu(self):
        """Initialize and configure the SMU"""
        # Reset and clear instrument
        self.smu.write('*RST')
        self.smu.write('TRIG:CLE')
        self.smu.write('TRACe:CLEAr')
        self.smu.write('*CLS')
        
        # Configure channels
        self.smu.write(':SYSTem:GROup:DEFine (1),(2)')
        
        # Configure output format
        self.smu.write(':FORMat:ELEMents:SENSe VOLT,CURRent,TIME')
        self.smu.write(':FORMat:DATA REAL,64')
        self.smu.write(':FORMat:BORDer SWAPped')
        
        # Configure output settings
        self.smu.write(':OUTPut1:ON:AUTO 1')
        self.smu.write(':OUTPut2:ON:AUTO 1')
        self.smu.write(':OUTPut1:FILTer:LPASs:STATe 1')
        self.smu.write(':OUTPut1:FILTer:LPASs:AUTO 1')
        self.smu.write(':OUTPut1:LOW GROund')
        self.smu.write(':OUTPut1:HCAPacitance:STATe 0')
        self.smu.write(':OUTPut2:FILTer:LPASs:STATe 1')
        self.smu.write(':OUTPut2:FILTer:LPASs:AUTO 1')
        self.smu.write(':OUTPut2:LOW GROund')
        self.smu.write(':OUTPut2:HCAPacitance:STATe 0')
        
        # Configure source settings
        self.smu.write(':SOURce1:FUNCtion:TRIGgered:CONTinuous 1')
        self.smu.write(':SOURce1:VOLTage:RANGe:AUTO 1')
        self.smu.write(':SOURce1:VOLTage:RANGe:AUTO:LLIMit MIN')
        self.smu.write(':SOURce1:VOLTage:LEVel:IMMediate:AMPLitude 0')
        self.smu.write(':SOURce1:VOLTage:LEVel:TRIGgered:AMPLitude 0')
        self.smu.write(':SOURce1:VOLTage:MODE FIXed')
        self.smu.write(':SENSe1:CURRent:DC:PROTection:LEVel 0.1')
        self.smu.write(':SOURce2:FUNCtion:MODE VOLTage')
        
        self.smu.write(':SOURce2:FUNCtion:TRIGgered:CONTinuous 1')
        self.smu.write(':SOURce2:VOLTage:RANGe:AUTO 1')
        self.smu.write(':SOURce1:VOLTage:RANGe:AUTO:LLIMit MIN')
        self.smu.write(':SOURce2:VOLTage:LEVel:IMMediate:AMPLitude 0')
        self.smu.write(':SOURce2:VOLTage:LEVel:TRIGgered:AMPLitude 0')
        self.smu.write(':SOURce2:VOLTage:MODE FIXed')
        self.smu.write(':SENSe2:CURRent:DC:PROTection:LEVel 0.1')
        self.smu.write(':SOURce2:FUNCtion:MODE VOLTage')
        
        # Configure digital I/O
        self.smu.write(':SOURce:DIGital:EXTernal14:FUNCtion DIO')
        self.smu.write(':SOURce:DIGital:DATA 0')
        
        # Configure trigger settings
        self.smu.write(':ARM1:TRANsient:LAYer:COUNt 1')
        self.smu.write(':TRIGger1:TRANsient:COUNt 1')
        self.smu.write(':ARM2:TRANsient:LAYer:COUNt 1')
        self.smu.write(':TRIGger2:TRANsient:COUNt 1')
        
        # Configure sense settings
        self.smu.write(':SENSe1:CURRent:DC:APERture {}'.format(self.current_measure_time))
        self.smu.write(':SENSe2:CURRent:DC:APERture {}'.format(self.current_measure_time))
        
        self.smu.write(':SENSe1:CURRent:DC:RANGe:AUTO ON')
        self.smu.write(':SENSe1:CURRent:DC:RANGe:AUTO:LLIMit MIN')
        self.smu.write(':SENSe1:CURRent:DC:RANGe:AUTO:MODE RES')
        self.smu.write(':SENSe1:CURRent:DC:RANGe:AUTO:THReshold MINimum')
        self.smu.write(':SENSe2:CURRent:DC:RANGe:AUTO ON')
        self.smu.write(':SENSe2:CURRent:DC:RANGe:AUTO:LLIMit MIN')
        self.smu.write(':SENSe2:CURRent:DC:RANGe:AUTO:MODE RES')
        self.smu.write(':SENSe2:CURRent:DC:RANGe:AUTO:THReshold MINimum')
        
        # Configure arm and trigger settings
        self.smu.write(':ARM1:TRANsient:LAYer:DELay 0.0')
        self.smu.write(':TRIGger1:TRANsient:DELay 0.05')
        self.smu.write(':TRIGger1:TRANsient:TIMer 0.1')
        self.smu.write(':TRIGger1:TRANsient:SOURce:SIGNal AINT')
        self.smu.write(':ARM2:TRANsient:LAYer:DELay 0.0')
        self.smu.write(':TRIGger2:TRANsient:DELay 0.05')
        self.smu.write(':TRIGger2:TRANsient:TIMer 0.1')
        self.smu.write(':TRIGger2:TRANsient:SOURce:SIGNal AINT')
        
        # Configure trigger output
        self.smu.write(':SOURce:DIGital:EXTernal9:FUNCtion DIO')
        self.smu.write(':SOURce:TOUTput:STATe 1')
        self.smu.write(':SOURce:TOUTput:SIGN EXT9')
        
        # Clear any errors
        self.smu.query(':SYSTem:ERRor:ALL?')
        self.smu.write('*CLS')
    
    @contextmanager
    def setup_plot(self, num, figsize=None):
        """Context manager for plot setup and display"""
        fig = plt.figure(num, figsize=figsize) if figsize else plt.figure(num)
        plt.ion()
        try:
            yield fig
        finally:
            plt.pause(0.01)
    
    def stabilize_spectrometer(self):
        """Run background stability test on spectrometer"""
        print('Running background stability test...')
        
        # Set integration time
        int_time_us = self.spec_integration_time * 1000
        self.spec.integration_time_micros(int_time_us)
        
        # Take initial reference spectra
        y_0_list = []
        for _ in range(self.spectra_repeats):
            x0 = self.spec.wavelengths()
            y0 = self.spec.intensities()
            y0 = convolve(y0, Box1DKernel(self.boxcar_width))
            y_0_list.append(y0)
        
        self.y0_data = np.mean(y_0_list, axis=0)
        self.x0_data = x0  # Save wavelengths
        
        # Check stability
        rms = float('inf')
        time_start = time.time()
        rms_array, time_array = [], []
        
        with self.setup_plot(1) as fig:
            while rms > self.optimal_rms and not self.interrupted:
                # Get new spectra
                noise_list = []
                for _ in range(self.spectra_repeats):
                    x1 = self.spec.wavelengths()
                    y1 = self.spec.intensities()
                    y1 = convolve(y1, Box1DKernel(self.boxcar_width))
                    noise_list.append(y1)
                
                noise = np.mean(noise_list, axis=0)
                array_noise = noise - self.y0_data
                
                # Calculate RMS
                peak_to_peak = np.amax(array_noise) - np.amin(array_noise)
                rms = peak_to_peak / 2 * np.sqrt(2)
                
                t = time.time() - time_start
                rms_array.append(rms)
                time_array.append(t)
                
                print(f'RMS: {rms:.2f}')
                
                # Plot noise
                plt.cla()
                plt.plot(x1, array_noise)
                plt.title(f'Spectrometer Noise - RMS: {rms:.2f}')
                plt.xlabel('Wavelength (nm)')
                plt.ylabel('Noise (counts)')
                
                # Update reference if not stabilized
                if rms > self.optimal_rms:
                    y_0_list = []
                    for _ in range(self.spectra_repeats):
                        x0 = self.spec.wavelengths()
                        y0 = self.spec.intensities()
                        y0 = convolve(y0, Box1DKernel(self.boxcar_width))
                        y_0_list.append(y0)
                    self.y0_data = np.mean(y_0_list, axis=0)
        
        print('Background is now stable.')
        
    def measure_sweep(self):
        """Perform voltage sweep and acquire measurements"""
        print('Starting measurement sweep...')
        
        # Prepare spectrometer
        spec_time_us = self.spec_integration_time * 1000
        self.spec.integration_time_micros(spec_time_us)
        
        # Reset SMU digital lines
        self.smu.write('*CLS')
        self.smu.write(':SOURce:DIGital:DATA 0')
        self.smu.write(':SOURce:DIGital:DATA 1')
        self.smu.write(':SOURce:DIGital:DATA 0')
        
        # Set initial voltage
        self.smu.write(':SOURce1:VOLTage:LEVel:IMMediate:AMPLitude 0')
        self.smu.write(':SOURce2:VOLTage:LEVel:IMMediate:AMPLitude 0')
        
        # Create voltage list for sweep
        voltage_list = np.arange(self.start_voltage, self.stop_voltage + self.step_size, self.step_size)
        
        # Setup measurement plot
        with self.setup_plot(1, (14, 10)) as fig:
            # Sweep through voltages
            for vg in voltage_list:
                if self.interrupted:
                    print("Measurement stopped by user")
                    break
                
                # Set voltage
                self.smu.write(f':SOURce1:VOLTage:LEVel:TRIG:AMPLitude {vg}')
                
                # Initialize and trigger measurement
                self.smu.write(':INIT (@1,2)')
                self.smu.write(':TRIGger:ALL (@1,2)')
                
                # Get measurement data
                ids = self.smu.query_binary_values(':SENSe2:DATA? CURRent,1', 'd', False)
                igs = self.smu.query_binary_values(':SENSe1:DATA? CURRent,1', 'd', False)
                
                # Get spectrum
                x_spec = self.spec.wavelengths()
                y_spec = self.spec.intensities()
                y_spec = convolve(y_spec, Box1DKernel(self.boxcar_width))
                y_data = y_spec - self.y0_data
                
                # Clear any errors
                self.smu.query(':SYSTem:ERRor:ALL?')
                self.smu.write('*CLS')
                
                # Print status
                print('=' * 40)
                print(f'Voltage: {vg:.2f} V')
                print(f'I_device: {igs[1]:.6e} A')
                print(f'J: {igs[1]/self.device_area:.6e} A/cm²')
                
                # Store data
                self.voltage.append(vg)
                self.current.append(igs[1])
                self.photocurrent.append(ids[1])
                self.current_density.append(igs[1]/self.device_area)
                
                # Process spectrum
                y_1nm = np.interp(self.cie_1nm_x, x_spec, y_data)
                int_time_norm = self.calibration_int_time_ms / self.spec_integration_time
                y_corrected = y_1nm * self.cal_y * int_time_norm
                
                self.spectra.append(y_data)
                self.interp_spectra.append(y_1nm)
                self.corrected_spectra.append(y_corrected)
                
                # Calculate EQE
                n_ph_nm = y_corrected * self.cie_1nm_x * 1e-9 / (constants.c * constants.h)
                eqe_val = np.sum(n_ph_nm) * 100 / (igs[1] * 6.24151E+18)
                
                if eqe_val < 0 or eqe_val > 100 or igs[1] < 0:
                    self.eqe.append(np.nan)
                else:
                    self.eqe.append(eqe_val)
                    print(f'EQE: {eqe_val:.2f} %')
                
                # Calculate luminance and efficiency
                cie_1nm_y_array = np.array(self.cie_1nm_y)
                luminosity_function = cie_1nm_y_array * 683  # lm/W
                lumens = y_corrected * luminosity_function
                self.lum_list_nm.append(lumens)
                
                lum = sum(lumens)
                self.luminance.append(lum)
                
                candelas = lum / np.pi
                brightness = candelas / self.device_area_m
                self.brightness.append(brightness)
                
                # Calculate power efficiency
                power_eff = lum / (vg * igs[1])
                if power_eff < 0 or igs[1] < 0:
                    self.power_efficiency.append(np.nan)
                else:
                    self.power_efficiency.append(power_eff)
                
                # Calculate current efficiency
                curr_eff = candelas / igs[1]
                if curr_eff < 0 or igs[1] < 0:
                    self.current_efficiency.append(np.nan)
                else:
                    self.current_efficiency.append(curr_eff)
                
                # Calculate color coordinates
                X = sum(y_corrected * self.cie_xyz_x)
                Y = sum(y_corrected * self.cie_xyz_y)
                Z = sum(y_corrected * self.cie_xyz_z)
                
                cie_x_val = X / (X + Y + Z)
                cie_y_val = Y / (X + Y + Z)
                
                self.cie_x.append(cie_x_val)
                self.cie_y.append(cie_y_val)
                
                # Update plots
                self.update_live_plots()
                
    def update_live_plots(self):
        """Update the real-time measurement plots"""
        # JV and Brightness
        ax1 = plt.subplot(321)
        ax1.cla()
        ax1.semilogy(self.voltage, self.current_density)
        ax1.set_xlabel('Voltage (V)')
        ax1.set_ylabel('Current density (A/cm²)')
        ax1.yaxis.label.set_color('blue')
        
        ax2 = ax1.twinx()
        ax2.semilogy(self.voltage, self.brightness, color='red')
        ax2.set_ylabel('Luminance (cd/m²)')
        ax2.yaxis.label.set_color('red')
        
        # EQE vs Brightness
        plt.subplot(323)
        plt.cla()
        plt.plot(self.brightness, self.eqe)
        plt.xlabel('Luminance (cd/m²)')
        plt.ylabel('EQE (%)')
        
        # JV curve
        plt.subplot(324)
        plt.cla()
        plt.plot(self.voltage, self.current_density)
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current density (A/cm²)')
        
        # Current spectrum
        plt.subplot(325)
        plt.cla()
        if self.spectra:
            x = self.spec.wavelengths()
            plt.plot(x, self.spectra[-1])
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Intensity (A.U.)')
        
        # EQE vs Voltage
        plt.subplot(122)
        plt.cla()
        plt.plot(self.voltage, self.eqe)
        plt.xlabel('Voltage (V)')
        plt.ylabel('EQE (%)')
        
    def calculate_mobility(self):
        """Calculate mobility and threshold voltage"""
        sqrt_ids = np.sqrt(np.abs(self.current))
        my_pwlf = pwlf.PiecewiseLinFit(self.voltage, sqrt_ids)
        breaks = my_pwlf.fit(3)
        
        # Find closest values
        abs_diff_func = lambda val: abs(val - breaks[2])
        closest_val = min(self.voltage, key=abs_diff_func)
        
        low_idx = self.voltage.index(closest_val)
        high_idx = self.voltage.index(min(self.voltage, key=lambda val: abs(val - breaks[3])))
        
        x_values = [closest_val, self.voltage[high_idx]]
        y_values = [sqrt_ids[low_idx], sqrt_ids[high_idx]]
        
        m, b = np.polyfit(x_values, y_values, 1)
        
        # Device parameters
        channel_length = 0.008  # cm
        channel_width = 0.3  # cm
        capacitance_pua = 15e-9  # F/cm²
        
        # Calculate threshold voltage and mobility
        v_threshold = -b / m
        mobility = (2 * channel_length) / (channel_width * capacitance_pua) * (m ** 2)
        
        return mobility, v_threshold
    
    def create_summary_plots(self):
        """Create summary plots of results"""
        # First summary plot
        with self.setup_plot(2) as fig:
            plt.subplot(321)
            plt.semilogy(self.voltage, self.current)
            plt.xlabel('Voltage (V)')
            plt.ylabel('log(I)')
            
            plt.subplot(323)
            plt.plot(self.voltage, self.eqe)
            plt.xlabel('Voltage (V)')
            plt.ylabel('EQE (%)')
            
            plt.subplot(325)
            plt.plot(self.cie_1nm_x, self.corrected_spectra[-1])
            plt.xlim([380, 750])
            plt.ylim(bottom=0)
            plt.ylabel('Intensity (A.U)')
            plt.xlabel('Wavelength (nm)')
            
            ax1 = plt.subplot(122)
            ax1.semilogy(self.voltage, self.current_density)
            ax1.set_xlabel('Voltage (V)')
            ax1.set_ylabel('Current density (A/cm²)')
            ax1.yaxis.label.set_color('blue')
            
            ax2 = ax1.twinx()
            ax2.semilogy(self.voltage, self.brightness, color='red')
            ax2.set_ylabel('Luminance (cd/m²)')
            ax2.yaxis.label.set_color('red')
        
        # Second summary plot
        with self.setup_plot(5, (20, 8)) as fig:
            plt.subplot(221)
            plt.plot(self.cie_1nm_x, self.corrected_spectra[-1])
            plt.xlim([380, 750])
            plt.ylim(bottom=0)
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Electro-Luminescence (A.U)')
            
            ax1 = plt.subplot(222)
            ax1.semilogy(self.voltage, self.current_density)
            ax1.set_xlabel('Voltage (V)')
            ax1.set_ylabel('Current density (A/cm²)')
            ax1.yaxis.label.set_color('blue')
            
            ax2 = ax1.twinx()
            ax2.semilogy(self.voltage, self.brightness, color='red')
            ax2.set_ylabel('Luminance (cd/m²)')
            ax2.yaxis.label.set_color('red')
            
            ax3 = plt.subplot(223)
            ax3.semilogy(self.current_density, self.eqe)
            ax3.set_xlabel('Current density (mA/cm²)')
            ax3.set_ylabel('EQE (%)')
            
            ax4 = plt.subplot(224)
            ax4.semilogy(self.brightness, self.eqe)
            ax4.set_xlabel('Luminance (cd/m²)')
            ax4.set_ylabel('EQE (%)')
            
    def save_data(self):
        """Save measurement data to CSV files"""
        print("Saving data to files...")
        
        # Save main results
        results_df = pd.DataFrame({
            'Voltage (V)': self.voltage,
            'Current (A)': self.current,
            'Current Density (A.cm^-2)': self.current_density,
            'Brightness (cd.m^-2)': self.brightness,
            'Current Efficiency': self.current_efficiency,
            'Power Efficiency': self.power_efficiency,
            'Photocurrent': self.photocurrent,
            'EQE (%)': self.eqe,
            'CIE_X': self.cie_x,
            'CIE_Y': self.cie_y
        })
        results_df.to_csv(self.filepath + self.filename_results, index=False)
        
        # Save wavelength data
        wavelength_df = pd.DataFrame({'Wavelength': self.cie_1nm_x})
        wavelength_df.to_csv(self.filepath + self.filename_wavelength, index=False)
        
        # Save spectral data
        spectra_df = pd.DataFrame(self.corrected_spectra)
        spectra_df.to_csv(self.filepath + self.filename_spectra, index=False)
        
        print(f"Data saved to {self.filepath}")
    
    def cleanup(self):
        """Clean up and close instruments"""
        try:
            # Reset and close SMU
            if hasattr(self, 'smu'):
                self.smu.write('*RST')
                self.smu.write('TRIG:CLE')
                self.smu.write('TRACe:CLEAr')
                self.smu.write('*CLS')
                self.smu.close()
            
            # Close spectrometer
            if hasattr(self, 'spec'):
                self.spec.close()
            
            # Close VISA resource manager
            if hasattr(self, 'rm'):
                self.rm.close()
            
            print("Instruments disconnected")
        except Exception as e:
            print(f"Error during cleanup: {e}")
    
    def run_measurement(self):
        """Execute the complete measurement sequence"""
        try:
            print("=" * 50)
            print("OLED CHARACTERIZATION SYSTEM")
            print("=" * 50)
            print("This program will measure OLED characteristics including:")
            print("- Current-voltage characteristics")
            print("- Luminance and spectral output")
            print("- Quantum efficiency")
            print("- Color coordinates")
            print("- Device mobility and threshold voltage")
            print("\nPlease ensure all instruments are connected properly.")
            print("=" * 50)
            
            # Load calibration data
            print("\nLoading calibration data...")
            self.load_calibration_data()
            
            # Connect to instruments
            print("\nConnecting to instruments...")
            self.connect_instruments()
            
            # Initialize SMU
            print("\nInitializing SMU...")
            self.initialize_smu()
            
            # Wait for user confirmation
            input("\nConnection established. Press Enter to continue...")
            
            # Stabilize spectrometer
            print("\nStabilizing spectrometer background...")
            self.stabilize_spectrometer()
            
            # Wait for user confirmation
            input("\nBackground is now stable. Press Enter to start experiment...")
            
            # Perform voltage sweep and measurements
            print("\nStarting voltage sweep measurements...")
            self.measure_sweep()
            
            # Calculate mobility and threshold voltage
            print("\nCalculating device parameters...")
            mobility, v_threshold = self.calculate_mobility()
            print(f"\nDevice Characteristics:")
            print(f"  - Mobility: {mobility:.5f} cm²/Vs")
            print(f"  - Threshold voltage: {v_threshold:.5f} V")
            
            # Create summary plots
            print("\nGenerating summary plots...")
            self.create_summary_plots()
            
            # Save data to files
            print("\nSaving measurement data...")
            self.save_data()
            
            # Display mobility and threshold results on final plot
            with self.setup_plot(7, (10, 6)) as fig:
                plt.clf()
                plt.text(0.5, 0.5, 
                       f"Device Characteristics:\n\n" +
                       f"Mobility: {mobility:.5f} cm²/Vs\n" +
                       f"Threshold Voltage: {v_threshold:.5f} V\n\n" +
                       f"Device Area: {self.device_area:.6f} cm²\n" +
                       f"Max Brightness: {max(self.brightness):.2f} cd/m²\n" +
                       f"Max EQE: {max(filter(lambda x: not np.isnan(x), self.eqe)):.2f} %",
                       ha='center', va='center', fontsize=14)
                plt.axis('off')
                plt.title("OLED Characterization Results", fontsize=16)
            
            # Block to keep plots visible until user closes them
            plt.ioff()
            plt.show(block=True)
            
        except Exception as e:
            print(f"\nError during measurement: {e}")
            traceback.print_exc()
        finally:
            # Always clean up resources
            self.cleanup()
            
            print("\n")
            print("         \\|||||/        ")
            print("         ( ^ ^ )         ")
            print("|--ooO-----(_)----------|")
            print("|                       |")
            print("|  OLED Characterization|")
            print("|      System v2.0      |")
            print("|------------------Ooo--|")
            print("         |__||__|        ")
            print("          ||  ||         ")
            print("         ooO  Ooo        ")