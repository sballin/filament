# filament

### To-do

Essential:

- Carry out fit_eigenmodes for one of the poloidal arrays

Cosmetic: 

- Combine signals, shotput, etc. into one file
- Replace VFR, OHR, etc with geometry.VFR, geometry.OHR...


### Code Structure

- __reconstruct_ip.py__
	- Calculates plasma current profile at each time point
	- Saves a movie of current profile evolution
- __fit_eigenmodes.py__
	- Calculates magnitude of eddy currents using eigenmodes
- __eigenmodes.py__
	- Calculates eigenmodes for stainless steel shells
- __geometry.py__
	- Establishes geometry of stainless steel shells
- __g_matrix.py__
	- Calculates Green's functions between all sensors and coils
	- Saves results in output/inductances(...).p to load in future
- __fields.pyx__
	- Can calculate magnetic field from a loop of current
	- Sped up using Cython, must run `make` before using 
- __shotput.py__
	- Something about loading shot data
- __signals.py__
	- Defines `Sensor` class and methods to obtain coil and sensor timeseries data
- __data_manipulation.py__
	- Methods for trimming, formatting timeseries data


### Eddy current magnitudes for unique sensors

![](resources/I_mags.png)


### Eigenmodes for a single shell

![](resources/eigenmodes.png)


### Initial current profile reconstruction

![]()

