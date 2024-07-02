import os
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import meep as mp
from meep import mpb
import pickle
import copy

class PhotonicCrystalWaveguideBands:
	c = 299792458
	field_components = {'x': 1, 'y': 0, 'z': 2}
	xz_component = None
	yz_component = None
	xy_component = None

	def __init__(self):
		pass

	def __str__(self):
		""" Return string specifying all relevant parameters
		"""
		return self.geometry_parameters['cross_section'] + '_bands_' + str(self.num_bands) + '_kstart_' + str(self.k_start) + '_kstop_' + str(self.k_stop) + '_kpoints_' + str(self.num_kpoints) + '_hx_' + str(self.geometry_parameters['hx']) + '_hy_' + str(self.geometry_parameters['hy']) + '_width_' + str(self.geometry_parameters['width']) + '_height_' + str(self.geometry_parameters['height']) + '_angle_' + str(self.geometry_parameters['angle']) + '_dboundary_' + str(self.geometry_parameters['dboundary']) + '_resolution_' + str(self.resolution)

	def setPath(self, path):
		""" Change working directory to path

		:param path: Directory to switch to
		"""
		self.path = path
		os.chdir(path)

	def buildGeometry(self, epsilon, hx, hy, width, dboundary, height=0.4, angle=np.pi/4, cross_section='rectangular', default_epsilon=1):
		"""Builds the goemetry, all dimensions are normalized to the unit cell spacing a

		:param epsilon: Material dielectric constant
		:param hx: Hole x (propagation direction) diameter
		:param hy: Hole y diameter
		:param width: Beam width
		:param height: Beam height (not used for triangular cross section, put something in anyway
		:param angle: Beam half-angle, only used for the triangular cross-section
		:param cross_section: 'rectangular' or 'triangular
		:param default_epsilon: Background material, set to air (epsilon = 1) by default
		"""

		self.geometry_lattice = mp.Lattice(size=mp.Vector3(1, 2, 2))
		discretization_points = 50
		theta = np.linspace(0, 2 * np.pi, discretization_points)
		vertices = [mp.Vector3(hx / 2 * np.cos(t), hy / 2 * np.sin(t), 0) for t in theta]
		if cross_section == 'rectangular':
			self.geometry = [mp.Block(center = mp.Vector3(0, 0, 0), size = mp.Vector3(2, width, height),
						  material = mp.Medium(epsilon = epsilon)),
					 mp.Prism(vertices, height = height + 2, center=mp.Vector3(0, 0, -height + 1),
						  material = mp.Medium(epsilon = default_epsilon))]
			self.geometry_lattice = mp.Lattice(size=mp.Vector3(1, width + 2 * dboundary, height + 2 * dboundary))
		elif cross_section == 'triangular':
			vertices_tri = [mp.Vector3(-1, -width / 2, 0.5), mp.Vector3(-1, width / 2, 0.5),
					mp.Vector3(-1, 0, 0.5 - width / (2 * np.tan(angle)))]
			self.geometry = [mp.Prism(vertices_tri, height = 2, axis = mp.Vector3(1, 0, 0),
						  material = mp.Medium(epsilon = epsilon)),
					 mp.Prism(vertices, height = width + 2, center=mp.Vector3(0, 0, -width + 1),
						  material = mp.Medium(epsilon = default_epsilon))]
			self.geometry_lattice = mp.Lattice(size=mp.Vector3(1, width + dboundary, width + dboundary))
		self.geometry_parameters = {'epsilon': epsilon, 'default_epsilon': default_epsilon, 'cross_section': cross_section,
					    'hx': hx, 'hy': hy, 'width': width, 'height': height, 'angle': angle, 'dboundary': dboundary}

	def setupSimulation(self, resolution, num_bands, k_start, k_stop, num_kpoints, default_epsilon=1):
		"""Sets simulation specific parameters

		:param resolution: Number of discretization points per axis
		:param num_bands: Number of Bands to calculate
		:param k_start: k-point to start simulation at
		:param k_stop: k-point to stop simulation at
		:param num_kpoints: Number of k-points to calculate the bandstructure for
		:param default_epsilon: Background material, set to air (epsilon = 1) by default
		"""

		self.resolution = resolution						# Spatial discretization points for geometric structure
		self.num_bands = num_bands						# Number of bands to solve for
		self.k_start = k_start							# First k-point
		self.k_stop = k_stop							# Last k-point
		self.num_kpoints = num_kpoints						# Number of k points to interpolate over
		k_points = [mp.Vector3(k_start, 0, 0), mp.Vector3(k_stop, 0, 0)]	# The k endpoints of interest for waveguide, this is 
											# usually the projected first Brillioun zone along
											# the propagation direction
		self.k_points = mp.interpolate(num_kpoints, k_points)
		self.te_ms = mpb.ModeSolver(num_bands=self.num_bands, k_points=self.k_points, geometry=self.geometry,
					    geometry_lattice=self.geometry_lattice, default_material=mp.Medium(epsilon = default_epsilon),
					    resolution=self.resolution)

		self.tm_ms = mpb.ModeSolver(num_bands=self.num_bands, k_points=self.k_points, geometry=self.geometry,
					    geometry_lattice=self.geometry_lattice, default_material=mp.Medium(epsilon = default_epsilon),
					    resolution=self.resolution)

	def run(self, save_plots=False, invert_band_diagram=False):
		"""Runs the simulation

		:param save_plots: A boolean specifying whether or not to save the standard plots.
		:return:
		"""
		self.fixed_field = 'E'
		self.unfixed_field = 'H'
		self.te_ms.run_yodd(mpb.fix_efield_phase)
		self.tm_ms.run_yeven(mpb.fix_efield_phase)
		self.te_freqs = self.te_ms.all_freqs
		self.tm_freqs = self.tm_ms.all_freqs
		self.num_kpoints = len(self.te_freqs)
		self.k = np.linspace(self.k_start, self.k_stop, self.num_kpoints)
		te_vg = self.te_ms.compute_group_velocities()
		tm_vg = self.tm_ms.compute_group_velocities()
		te_vg = np.array([[i + 1, te_vg[i].x, te_vg[i].y, te_vg[i].y] for i in range(len(te_vg))])
		tm_vg = np.array([[i + 1, tm_vg[i].x, tm_vg[i].y, tm_vg[i].y] for i in range(len(tm_vg))])
		columns = ['Band Number', '$v_{gx}$', '$v_{gy}$', '$v_{gz}$']
		self.te_vg_df = pd.DataFrame(te_vg, columns=columns)
		self.tm_vg_df = pd.DataFrame(tm_vg, columns=columns)
		self.epsilon_geometry = self.tm_ms.get_epsilon()
		self.ny, self.nx, self.nz = self.epsilon_geometry.shape
		#self.te_ms.k_points = [mp.Vector3(0.5, 0, 0)]
		#self.te_ms.k_points = [mp.Vector3(0.5, 0, 0)]
		#self.te_ms.run_te(mpb.fix_efield_phase)
		#self.tm_ms.run_tm(mpb.fix_efield_phase)
		H_te = []
		D_te = []
		E_te = []
		S_te = []
		H_tm = []
		D_tm = []
		E_tm = []
		S_tm = []
		for n in range(self.num_bands):
			H_te.append(self.te_ms.get_hfield(n+1, bloch_phase=True))
			D_te.append(self.te_ms.get_dfield(n+1, bloch_phase=True))
			E_te.append(self.te_ms.get_efield(n+1, bloch_phase=True))
			S_te.append(self.te_ms.get_poynting(n+1))
			H_tm.append(self.tm_ms.get_hfield(n+1, bloch_phase=True))
			D_tm.append(self.tm_ms.get_dfield(n+1, bloch_phase=True))
			E_tm.append(self.tm_ms.get_efield(n+1, bloch_phase=True))
			S_tm.append(self.tm_ms.get_poynting(n+1))
		self.fields_data_te = {'H': H_te, 'D': D_te, 'E': E_te, 'S': S_te}
		self.fields_data_tm = {'H': H_tm, 'D': D_tm, 'E': E_tm, 'S': S_tm}
		self.calcBandGaps(condition='X')
		self.getEnergyDistribution()
		if save_plots:
			self.plotBands(polarization='all', gaps='one', units=None, invert=False, save=True)
			self.getField('Hz', 'TE', 1)
			fig1 = self.plotSlices(0, 0, 0)
			plt.savefig('Hz_TE_band1_' + str(self) + '.png')
			self.getField('Hz', 'TE', 2)
			fig2 = self.plotSlices(0, 0, 0)
			plt.savefig('Hz_TE_band2_' + str(self) + '.png')

	def saveData(self):
		"""Saves the simulated data to data.p in path
		"""
		data = {'geometry_parameters': self.geometry_parameters,
			'resolution': self.resolution,
			'num_bands': self.num_bands,
			'num_kpoints': self.num_kpoints,
			'k': self.k,			
			'te_freqs': self.te_freqs,
			'tm_freqs': self.tm_freqs,
			'te_vg_df': self.te_vg_df,
			'tm_vg_df': self.tm_vg_df,
			'te_quasi_gaps_df': self.te_quasi_gaps_df,
			'tm_quasi_gaps_df': self.tm_quasi_gaps_df,
			'epsilon_geometry': self.epsilon_geometry,
			'fields_data_te': self.fields_data_te,
			'fields_data_tm': self.fields_data_tm,
			'field_energy': self.field_energy,
			'power_fractions_df': self.power_fractions_df}
		with open('data.p', 'wb') as f:
			pickle.dump(data, f)

	def loadData(self):
		"""Loads already simulated data from data.p in path and recreates all data structures with the exception of the
		modesolvers not being able to calculate anything
		"""
		data = pickle.load(open('data.p','rb'))
		self.geometry_parameters = data['geometry_parameters']
		self.buildGeometry(self.geometry_parameters['epsilon'], 
				   self.geometry_parameters['hx'],
				   self.geometry_parameters['hy'],
				   self.geometry_parameters['width'],
				   self.geometry_parameters['dboundary'],
				   height=self.geometry_parameters['height'],
				   angle=self.geometry_parameters['angle'],
				   cross_section=self.geometry_parameters['cross_section'],
				   default_epsilon=self.geometry_parameters['default_epsilon'])
		self.resolution = data['resolution']
		self.num_bands = data['num_bands']
		self.num_kpoints = data['num_kpoints']
		self.k = data['k']
		self.k_start = self.k[0]
		self.k_stop = self.k[-1]
		self.setupSimulation(self.resolution, self.num_bands, self.k_start, self.k_stop, self.num_kpoints)
		self.te_freqs = data['te_freqs']
		self.tm_freqs = data['tm_freqs']
		self.te_vg_df = data['te_vg_df']
		self.tm_vg_df = data['tm_vg_df']
		self.te_quasi_gaps_df = data['te_quasi_gaps_df']
		self.tm_quasi_gaps_df = data['tm_quasi_gaps_df']
		self.epsilon_geometry = data['epsilon_geometry']
		self.fields_data_te = data['fields_data_te']
		self.fields_data_tm = data['fields_data_tm']
		self.field_energy = data['field_energy']
		self.power_fractions_df = data['power_fractions_df']
		self.ny, self.nx, self.nz = self.epsilon_geometry.shape
		self.calcBandGaps(condition='X')
		self.fixed_field = 'E'
		self.unfixed_field = 'H'

	def calcBandGaps(self, condition='X'):
		"""Calculates all quasi-band gaps (considering only points below the light line) present in the bandstructure

		:param condition: 'X' for considering only the X point, 'all' for considering all points below the light line,
		if interested in full band gaps including points above the light line, this is output directly by the simulation
		"""
		te_valid_band_range_data = np.zeros((self.num_bands, 2))
		tm_valid_band_range_data = np.zeros((self.num_bands, 2))
		te_quasi_gaps = []
		tm_quasi_gaps = []

		for n, te_freq, tm_freq in zip(np.arange(0, self.num_bands), np.transpose(self.te_freqs), np.transpose(self.tm_freqs)):
			# Find all points below the light line
			te_valid_freqs = []
			tm_valid_freqs = []
			for i in range(len(self.k)):
				if te_freq[i] < self.k[i]:
					te_valid_freqs.append(te_freq[i])
				if tm_freq[i] < self.k[i]:
					tm_valid_freqs.append(tm_freq[i])

			if condition == 'X':
				# Consider only the X point
				if len(te_valid_band_range_data) != 0 and len(te_valid_freqs) != 0:
					te_valid_band_range_data[n,0] = te_valid_freqs[-1]
					te_valid_band_range_data[n,1] = te_valid_freqs[-1]
				if len(tm_valid_band_range_data) != 0 and len(tm_valid_freqs) != 0:
					tm_valid_band_range_data[n,0] = tm_valid_freqs[-1]
					tm_valid_band_range_data[n,1] = tm_valid_freqs[-1]
			elif condition == 'all':
				# Consider all points below the light line
				if len(te_valid_band_range_data) != 0 and len(te_valid_freqs) != 0:
					te_valid_band_range_data[n,0] = np.min(te_valid_freqs)
					te_valid_band_range_data[n,1] = np.max(te_valid_freqs)
				if len(tm_valid_band_range_data) != 0 and len(tm_valid_freqs) != 0:
					tm_valid_band_range_data[n,0] = np.min(tm_valid_freqs)
					tm_valid_band_range_data[n,1] = np.max(tm_valid_freqs)

		# Calculate relevant band gap data (midgap, gap size, upper/lower band and band edge), neglect bands smaller
		# than 1% as they are likely numerical error
		for n in range(self.num_bands - 1):
			te_midgap = (te_valid_band_range_data[n+1,0] + te_valid_band_range_data[n,1]) / 2
			tm_midgap = (tm_valid_band_range_data[n+1,0] + tm_valid_band_range_data[n,1]) / 2
			te_gap_size = 100 * (te_valid_band_range_data[n+1,0] - te_valid_band_range_data[n,1]) / te_midgap
			tm_gap_size = 100 * (tm_valid_band_range_data[n+1,0] - tm_valid_band_range_data[n,1]) / tm_midgap
			if te_gap_size > 1:
				te_quasi_gaps.append([te_gap_size, n + 1, te_valid_band_range_data[n,1], n + 2,
						      te_valid_band_range_data[n+1,0], te_midgap])
			if tm_gap_size > 1:
				tm_quasi_gaps.append([tm_gap_size, n + 1, tm_valid_band_range_data[n,1], n + 2,
						      tm_valid_band_range_data[n+1,0], tm_midgap])
			columns=['Band Gap Size (%)', 'Lower Band Number', 'Lower Band Edge ($\omega$a/2$\pi$c)',
				 'Upper Band Number', 'Upper Band Edge ($\omega$a/2$\pi$c)', 'Midgap ($\omega$a/2$\pi$c)']
			self.te_quasi_gaps = te_quasi_gaps
			self.tm_quasi_gaps = tm_quasi_gaps
			if len(self.te_quasi_gaps) != 0:
				self.te_quasi_gaps_df = pd.DataFrame(np.array(te_quasi_gaps), columns=columns)
			else:
				self.te_quasi_gaps_df = []
			if len(self.tm_quasi_gaps) != 0:
				self.tm_quasi_gaps_df = pd.DataFrame(np.array(tm_quasi_gaps), columns=columns)
			else:
				self.te_quasi_gaps_df = []

	def plotBands(self, polarization='all', gaps='one', a=None, units = None, invert=False, save=False):
		"""Plots the bandstructure

		:param polarization: 'TE', 'TM', or 'all'
		:param gaps: 'one' plots just the first band gap, 'all' plots them all
		:param a: Unit cell length in nm
		:param units: None for unitless, or 'THz'
		:param invert: Inverts the x (k) axis if True
		:param save: Save the plot
		"""
		if units is None:
			nu_constant = 1
			k_constant = 1
			k = self.k
			xlabel = 'Wave vector (ka/2$\pi$)'
			ylabel = 'Frequency (' + r'$\nu$' + 'a/c)'
		elif units == 'THz':
			if a is None:
				raise ValueError('a must be specified')
			nu_constant = 1e-3 * self.c / a
			k_constant = 1000 * 2 * np.pi / a
			xlabel = 'Wave vector k (rad/$\mu$m)'
			ylabel = 'Frequency ' +r'$\nu$' + ' (THz)'

		te_freqs = nu_constant * np.transpose(self.te_freqs)
		tm_freqs = nu_constant * np.transpose(self.tm_freqs)		
		k = k_constant * self.k
		fig, ax = plt.subplots()
		custom_lines = []
		leg = []
		if polarization == 'TE' or polarization == 'all':
			min_te_freq = np.min(te_freqs)
			max_te_freq = np.max(te_freqs)
			min_freq = min_te_freq
			max_freq = max_te_freq
			custom_lines.append(Line2D([0], [0], color='red', lw=2))
			leg.append('TE Bands')
			for te_freq in te_freqs:
				ax.scatter(k, te_freq, color='red', facecolors = 'none')
				ax.plot(k, te_freq, color='red')
			if gaps == 'one' and len(self.te_quasi_gaps) != 0:
				te_gap = self.te_quasi_gaps[0]
				ax.fill_between(k, nu_constant * te_gap[2], nu_constant * te_gap[4], color='red', alpha=0.2)
			elif gaps == 'all' and len(self.te_quasi_gaps) != 0:
				for te_gap in self.te_quasi_gaps:
					ax.fill_between(k, nu_constant * te_gap[2], nu_constant * te_gap[4], color='red', alpha=0.2)
			else:
				print('No TE Band Gaps or invalid gaps choice')
		if polarization == 'TM' or polarization == 'all':
			min_tm_freq = np.min(tm_freqs)
			max_tm_freq = np.max(tm_freqs)
			min_freq = min_tm_freq
			max_freq = max_tm_freq
			custom_lines.append(Line2D([0], [0], color='blue', lw=2))
			leg.append('TM Bands')
			for tm_freq in tm_freqs:
				ax.scatter(k, tm_freq, color='blue')
				ax.plot(k, tm_freq, color='blue')
			if gaps == 'one' and len(self.tm_quasi_gaps) != 0:
				tm_gap = self.tm_quasi_gaps[0]
				ax.fill_between(k, nu_constant * tm_gap[2], nu_constant * tm_gap[4], color='blue', alpha=0.2)
			elif gaps == 'all' and len(self.tm_quasi_gaps) != 0:
				for tm_gap in self.tm_quasi_gaps:
					tm_gap = nu_const * tm_gap
					ax.fill_between(k, nu_constant * tm_gap[2], nu_constant * tm_gap[4], color='blue', alpha=0.2)
			else:
				raise ValueError('No TM Band Gaps or Invalid gaps choice')
		if polarization == 'all':
			min_freq = np.min([min_te_freq, min_tm_freq])
			max_freq = np.max([max_te_freq, max_tm_freq])
		ax.set_ylim([min_freq, max_freq])
		ax.set_xlim([k[0], k[-1]])		
		ax.set_xlabel(xlabel, size=16)
		ax.set_ylabel(ylabel, size=16)
		ax.grid(True)
		ax.fill_between(k, nu_constant / k_constant * k, max_freq, color = 'purple', alpha= 0.7)
		ax.legend(custom_lines, leg, loc='lower right')
		if invert:
			ax.invert_xaxis()
		if save:
			plt.savefig('bandstructure_' + str(self) + '.png')

	def listGaps(self, polarization = 'all', gaps='one', a=None, units=None):
		"""Lists band gaps
 
		:param polarization: 'TE', 'TM', or 'all'
		:param gaps: 'one' plots just the first band gap, 'all' plots them all
		:param a: Unit cell length in nm
		:param units: None for unitless, 'nm', or 'THz'
		"""
		te_quasi_gaps_df_temp = copy.deepcopy(self.te_quasi_gaps_df)
		tm_quasi_gaps_df_temp = copy.deepcopy(self.tm_quasi_gaps_df)
		if units is None:
			pass
		elif units == 'THz' and a is not None:
			constant = 1e-3 * self.c / a
			new_columns = ['Band Gap Size (%)', 'Lower Band Number', 'Lower Band Edge (THz)',
				       'Upper Band Number', 'Upper Band Edge (THz)', 'Midgap (THz)']
			te_quasi_gaps_df_temp.columns = new_columns
			tm_quasi_gaps_df_temp.columns = new_columns
			te_quasi_gaps_df_temp['Lower Band Edge (THz)'] = te_quasi_gaps_df_temp['Lower Band Edge (THz)'] * constant
			te_quasi_gaps_df_temp['Upper Band Edge (THz)'] = te_quasi_gaps_df_temp['Upper Band Edge (THz)'] * constant
			te_quasi_gaps_df_temp['Midgap (THz)'] = te_quasi_gaps_df_temp['Midgap (THz)'] * constant
			tm_quasi_gaps_df_temp['Lower Band Edge (THz)'] = tm_quasi_gaps_df_temp['Lower Band Edge (THz)'] * constant
			tm_quasi_gaps_df_temp['Upper Band Edge (THz)'] = tm_quasi_gaps_df_temp['Upper Band Edge (THz)'] * constant
			tm_quasi_gaps_df_temp['Midgap (THz)'] = tm_quasi_gaps_df_temp['Midgap (THz)'] * constant
		elif units == 'nm' and a is not None:
			new_columns = ['Band Gap Size (%)', 'Lower Band Number', 'Lower Band Edge (nm)',
				       'Upper Band Number', 'Upper Band Edge (nm)', 'Midgap (nm)']
			te_quasi_gaps_df_temp.columns = new_columns
			tm_quasi_gaps_df_temp.columns = new_columns
			te_quasi_gaps_df_temp['Lower Band Edge (nm)'] = a / te_quasi_gaps_df_temp['Lower Band Edge (nm)']
			te_quasi_gaps_df_temp['Upper Band Edge (nm)'] = a / te_quasi_gaps_df_temp['Upper Band Edge (nm)']
			te_quasi_gaps_df_temp['Midgap (nm)'] = a / te_quasi_gaps_df_temp['Midgap (nm)']
			tm_quasi_gaps_df_temp['Lower Band Edge (nm)'] = a / tm_quasi_gaps_df_temp['Lower Band Edge (nm)']
			tm_quasi_gaps_df_temp['Upper Band Edge (nm)'] = a / tm_quasi_gaps_df_temp['Upper Band Edge (nm)']
			tm_quasi_gaps_df_temp['Midgap (nm)'] = a / tm_quasi_gaps_df_temp['Midgap (nm)']
		elif (units == 'THz' or units == 'nm') and a is not None:
			raise ValueError('a must be specified')
		else:
			raise ValueError('Invalid units')
		if len(self.te_quasi_gaps) != 0 and (polarization == 'TE' or polarization == 'all'):
			if gaps == 'one':
				print('TE Band Gap')
				display(te_quasi_gaps_df_temp.loc[[0]])
			elif gaps == 'all':
				print('TE Band Gaps')
				display(te_quasi_gaps_df_temp)
		else:
			print('No TE Band Gaps')
		if len(self.tm_quasi_gaps) != 0 and (polarization == 'TM' or polarization == 'all'):
			if gaps == 'one':
				print('TM Band Gap')
				display(tm_quasi_gaps_df_temp.loc[[0]])
			if gaps == 'all':
				print('TM Band Gaps')
				display(tm_quasi_gaps_df_temp)

	def getEnergyDistribution(self):
		"""Calculates the energy percentage in the real and imaginary part of each field component
		"""
		field_energy = []
		for band in range(self.num_bands):		
			self.te_ms.get_hfield(band + 1, bloch_phase=True)
			field_energy.append(self.te_ms.compute_field_energy()[1:])
		columns = ['Re($' + self.unfixed_field + '_x$) %', 'Im($' + self.unfixed_field + '_x$) %',
			   'Re($' + self.unfixed_field + '_y$) %', 'Im($' + self.unfixed_field + '_y$) %',
			   'real($' + self.unfixed_field + '_z$) %', 'Im($' + self.unfixed_field + '_z$) %']
		self.field_energy = np.array(field_energy)
		self.power_fractions_df = pd.DataFrame(100 * self.field_energy, columns=columns)
		self.power_fractions_df.insert(0,'Band Number', np.arange(0, self.num_bands) + 1, True)

	def getField(self, field, polarization, band, part='max'):
		"""Gets the mode field profiles at the X point

		:param field: Allowed vector fields are 'H', 'D', 'E', and 'S'
		Scalar fields specified by appending the component (x, y, z, m - magnitide),
		e.g. 'Hx', 'Hy', 'Hz', 'Hm'
		Energy desnity of the component can be found by then appending E e.g. 'HxE', 'HyE', 'HzE', 'HmE'
		:param polarization: Allowed polarizations are 'TE' and 'TM'
		:param band: Allow bands range from 1 to num_bands
		:param part: 'real', 'imaginary', 'magnitude', 'argument', or 'max'
		Unfortunatuely the field component phase is arbitrary (though not random), for the sake of plotting,
		for a given band these simulations fix either E or H to be completely real meaninging the other will
		be purely imaginary. 'max' chooses whichever of the real of imaginary part is larger
		"""
		self.field = field
		if polarization == 'TE':
			field_data = self.fields_data_te[field[0]][band - 1]
		elif polarization == 'TM':
			field_data = self.fields_data_tm[field[0]][band - 1]
		else:
			raise ValueError("Invalid Polarization, must be 'TE' or 'TM'")

		if part == 'real':
			field_data = np.real(field_data)
		elif part == 'imaginary':
			field_data = np.imag(field_data)
		elif part == 'magnitude':
			field_data = np.abs(field_data)
		elif part == 'argument':
			field_data = np.angle(field_data)
		elif part == 'max':
			temp = np.reshape(np.array(self.field_energy[band - 1]),(3,2))
			is_imaginary = np.where(temp == np.amax(temp))[1] == 1
			if (is_imaginary and field[0] != self.fixed_field) or (not is_imaginary and field[0] == self.fixed_field):
				field_data = np.imag(field_data)
			elif (not is_imaginary and field[0] != self.fixed_field) or (is_imaginary and field[0] == self.fixed_field):
				field_data = np.real(field_data)
			else:
				raise ValueError('Huh')
		else:
			raise ValueError('Invalid part')

		try:
			self.field_type = 'scalar'
			if field[1] == 'x':
				field_data = field_data[...,1].reshape(self.ny, self.nx, self.nz)
			elif field[1] == 'y':
				field_data = field_data[...,0].reshape(self.ny, self.nx, self.nz)
			elif field[1] == 'z':
				field_data = field_data[...,2].reshape(self.ny, self.nx, self.nz)
			elif field[1] == 'm':
				field_data = np.linalg.norm(field_data, axis = -1).reshape(self.ny, self.nx, self.nz)
			else:
				raise ValueError("Invalid Field Component, must be 'x' or 'y', 'z', or 'm'")
			try:
				if field[2] == 'E':
					if field[0] == 'H':
						field_data = field_data ** 2
					elif field[0] == 'D':
						field_data = field_data ** 2 / self.epsilon_geometry
					elif field[0] == 'E':
						field_data = field_data ** 2 * self.epsilon_geometry
					elif field[0] == 'S':
						pass
			except:
				pass
		except IndexError:
			self.field_type = 'vector'
			field_data = field_data.reshape(self.ny, self.nx, self.nz, 3)

		self.field_data = field_data

	def plotSlices(self, x, y, z):
		"""Places 2D slices of the current field_data in the xz, yz, and xy planes

		:param x: yz plan cut coordinate (0 is the center of the simulation)
		:param y: xz plan cut coordinate (0 is the center of the simulation)
		:param z: xy plan cut coordinate (0 is the center of the simulation)
		"""
		if len(self.field_data.shape) == 3:
			field_xz = self.field
			field_yz = self.field
			field_xy = self.field
			field_data_xz = self.field_data
			field_data_yz = self.field_data
			field_data_xy = self.field_data
		elif len(self.field_data.shape) == 4 and self.xz_component != None and self.yz_component != None and self.xy_component != None:
			field_xz = self.field + '$_' + self.xz_component + '$'
			field_yz = self.field + '$_' + self.yz_component + '$'
			field_xy = self.field + '$_' + self.xy_component + '$'
			field_data_xz = self.field_data[:, :, :, self.field_components[self.xz_component]]
			field_data_yz = self.field_data[:, :, :, self.field_components[self.yz_component]]
			field_data_xy = self.field_data[:, :, :, self.field_components[self.xy_component]]
		elif len(self.field_data.shape) == 4 and (self.xz_component == None or self.yz_component == None or self.xy_component == None):
			raise ValueError('Slice field components must be specified')
		else:
			raise ValueError('Invalid field')

		fig = plt.figure()
		plt.subplot(131)
		xz_max = np.max(np.abs(field_data_xz))
		plt.title(field_xz + ' xz slice')
		plt.contour(self.epsilon_geometry[:, x + self.nx//2, :].T, cmap='binary')
		plt.imshow(field_data_xz[:, x + self.nx//2, :].T, interpolation='spline36', cmap='RdBu_r', alpha=0.9,
			   vmin = -xz_max, vmax = xz_max)
		ax = plt.gca()
		ax.invert_yaxis()
		plt.axis('off')
		plt.colorbar(orientation = 'horizontal')

		plt.subplot(132)
		yz_max = np.max(np.abs(field_data_yz))
		plt.title(field_yz + ' yz slice')
		plt.contour(self.epsilon_geometry[y + self.ny//2, :, :].T, cmap='binary')
		plt.imshow(field_data_yz[y + self.ny//2, :, :].T, interpolation='spline36', cmap='RdBu_r', alpha=0.9,
			   vmin = -yz_max, vmax = yz_max)
		ax = plt.gca()
		ax.invert_yaxis()
		plt.axis('off')
		plt.colorbar(orientation = 'horizontal')

		plt.subplot(133)
		xy_max = np.max(np.abs(field_data_xy))
		plt.title(field_xy + ' xy slice')
		plt.contour(self.epsilon_geometry[:, :, -z + self.nz//2].T, cmap='binary')
		plt.imshow(field_data_xy[:, :, z + self.nz//2].T, interpolation='spline36', cmap='RdBu_r', alpha=0.9,
			   vmin = -xy_max, vmax = xy_max)
		ax = plt.gca()
		ax.invert_yaxis()
		plt.axis('off')
		plt.colorbar(orientation = 'horizontal')

