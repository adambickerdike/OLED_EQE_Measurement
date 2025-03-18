# OLED_EQE_Measurement

# OLED Characterisation System

## Overview
This Python-based system facilitates comprehensive characterisation of Organic Light-Emitting Diodes (OLEDs) using SCPI command interfaces. It provides automated measurement of current-voltage characteristics, luminance, spectral output, efficiency metrics, and colour coordinates.

## Features
- Automated voltage sweep measurements
- Real-time data visualisation
- Calculation of key device parameters:
  - Current density
  - Luminance
  - External Quantum Efficiency (EQE)
  - Current and power efficiency
  - CIE colour coordinates
  - Charge carrier mobility
  - Threshold voltage
- Data export to CSV format
- Comprehensive graphical analysis

## Requirements

### Hardware
- Keysight B2912A Precision SMU (or compatible SMU with SCPI support)
- Ocean Optics/Seabreeze-compatible spectrometer
- GPIB interface adaptor
- Integrating sphere setup

### Software Dependencies
- Python 3.6 or higher
- NumPy
- Pandas
- Matplotlib
- SciPy
- Astropy
- PyVISA
- PWLF (Piecewise Linear Fitting)
- Seabreeze

Install dependencies using:
```bash
pip install numpy pandas matplotlib scipy astropy pyvisa pwlf seabreeze
```

### Calibration Files
Ensure the following calibration files are available in the same directory:
- `CIE_1nm_NIR_test.csv` - CIE colour matching functions
- `1000um_Calibration_Curve_27_03_2022.csv` - Spectrometer calibration data
- `CIE_xyz.csv` - CIE tristimulus functions

## Usage

### Configuration
Before running a measurement, adjust the device parameters in the `__init__` method:

```python
# Device parameters
self.dev_width = 0.025  # cm - update with your device's width
self.dev_length = 0.01  # cm - update with your device's length

# Measurement settings
self.start_voltage = 0.0  # V - starting voltage for sweep
self.stop_voltage = 10.0  # V - ending voltage for sweep
self.step_size = 0.5     # V - voltage step size

# Output file paths
self.filepath = 'D:\\Integrating_Sphere_Results\\06_09_2024\\'  # your directory
self.filename_results = 'OLED_Results.csv'  # your result filename
```

### Running a Measurement
1. Connect your OLED device to the SMU (typically channel 1 for voltage application)
2. Place the device in the integrating sphere with the spectrometer
3. Run the script from your IDE or command line:
   ```
   python oled_measurement.py
   ```
4. Follow the on-screen prompts:
   - First, the script will establish connections to instruments
   - It will then stabilise the spectrometer background
   - Finally, it will perform the voltage sweep and measurements

### Output
The system generates three CSV files:
1. Main results file: Contains voltage, current, brightness, efficiency metrics, and colour coordinates
2. Wavelength file: Contains the wavelength values for spectral data
3. Spectral data file: Contains the corrected spectra for each measurement point

Additionally, the system provides several visualisation plots:
- Current-voltage characteristics with brightness overlay
- EQE vs. brightness
- EQE vs. voltage
- Emission spectrum
- Device parameters summary

## Troubleshooting

### Common Issues
- **Connection errors**: Ensure GPIB addresses match and all cables are secure
- **Calibration file errors**: Verify file paths and proper CSV format
- **Spectrometer timeout**: Reduce integration time if the spectrometer isn't responding
- **Poor EQE values**: Check for ambient light contamination or improper device placement

### Error Handling
The system includes robust error handling to prevent catastrophic failures. If an error occurs during measurement, the system will:
1. Display an error message
2. Safely close all instrument connections
3. Save any collected data where possible

## Customisation
The system can be extended for different measurement types:
- For transient measurements, modify the trigger settings in `initialize_smu()`
- For temperature-dependent measurements, add an external temperature controller interface
- For automated device testing, implement a device handler class

## Acknowledgements
by adam b
