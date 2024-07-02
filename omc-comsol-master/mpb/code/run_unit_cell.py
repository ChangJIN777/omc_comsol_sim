from photonic_crystal_waveguide import PhotonicCrystalWaveguideBands as PC
import numpy as np
import os

# Path to operate in and save files to
base_path = '/mnt/12B92FAC0F9C831C/Simulations_and_Scripts/mpb/data/'  # Very important to have a slash at the end of this path

# Physical Parameters:
# Note that all parameters are normalized to the unit cell length 'a', this is fundamental to how this code is written and possibly how this mode solver operates
epsilon = 5.779						# Material dielectric constant (n**2 = 2.404**2 = 5.779 for diamond)
default_epsilon = 1					# Background material dielectric constant12
cross_section = 'triangular'				# Beam cross-section
hx_sweep = np.linspace(0.431, 0.4, 1)			# Hole y diameter
hy_sweep = np.linspace(1.017, 0.4, 1)			# Hole x diameter
width_sweep = np.linspace(1.602, 2, 1)			# Beam width
height_sweep = np.linspace(0.505, 0.5, 1)		# Beam height, not used for triangular cross-section simulations
angle_sweep = np.linspace(0.611, 0.9, 1)		# Beam half-angle in radians, not used for rectangular cross-section simulations
dboundary_sweep = np.linspace(10, 5, 1)			# Distance to simulation periodic boundary, very important for points above the light line, not important
							# at the X point (for the lower bands of interest negligible ~3e-4 error at 1, fully converged to within
							# numerical error at 3)

# Simulation Parameters
resolution_sweep = np.linspace(10, 20, 1)		# Unit cell discretization points (just of dielectric material), RMS error is ~1e-4 at resolution = 10
							# and decreases at ~10 dB/ 7 points, marginal gains by increasing this
num_bands = 4						# Number of bands to calculate
k_start = 0.3						# Starting k point (unitless)
k_stop = 0.5						# Ending k points (max 0.5, X point)
num_kpoints = 20					# Number of k points

for hx in list(hx_sweep):
	for hy in list(hy_sweep):
		for width in list(width_sweep):
			for dboundary in list(dboundary_sweep):
				for height in list(height_sweep):
					for angle in (angle_sweep):
						for resolution in list(resolution_sweep):
							os.chdir(base_path)
							pc = PC()
							pc.buildGeometry(epsilon, float(hx), float(hy), float(width), float(dboundary), height=float(height), angle=float(angle), cross_section=cross_section, default_epsilon=default_epsilon)
							pc.setupSimulation(resolution, num_bands, k_start, k_stop, num_kpoints, default_epsilon=default_epsilon)
							folder_name = str(pc)
							os.mkdir(folder_name)
							pc.setPath(base_path + folder_name)
							pc.run(save_plots=True, invert_band_diagram=False)
							pc.saveData()
							del pc
