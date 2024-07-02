from photonic_crystal_waveguide import PhotonicCrystalWaveguideBands as PC
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)

class Sweep:
	columns = ['bands', 'kstart', 'kstop', 'kpoints', 'hx', 'hy', 'width', 'height', 'angle', 'dboundary', 'resolution']
	def __init__(self, directory):
		self.directory = directory
		parameters = []
		self.te_band_freqs = []
		self.tm_band_freqs = []
		self.te_band_gaps = []
		self.tm_band_gaps = []
		for folder in sorted_alphanumeric(os.listdir(directory)):
			parameters.append([float(n) for n in folder.split('_')[::2][1:]])
			os.chdir(directory + folder)
			pc = PC()
			pc.loadData()
			self.te_band_freqs.append(pc.te_freqs)
			self.tm_band_freqs.append(pc.tm_freqs)
			self.te_band_gaps.append(pc.te_quasi_gaps_df)
			self.tm_band_gaps.append(pc.tm_quasi_gaps_df)

		self.parameters = pd.DataFrame(parameters, columns=self.columns)
		self.num_kpoints = pc.num_kpoints
		self.k = pc.k
		self.num_bands = pc.num_bands
		self.bands = np.arange(1, self.num_bands + 1)


	def plotConvergence(self, k_ind, legend, log):
		K, DBOUNDARY = np.meshgrid(self.k, np.linspace(1,10,9))
		te_relative_error = [np.abs(self.te_band_freqs[i+1] - self.te_band_freqs[i]) for i in range(len(self.te_band_freqs)-1)]
		te_absolute_error = [np.abs(self.te_band_freqs[i] - self.te_band_freqs[-1]) for i in range(len(self.te_band_freqs)-1)]
		tm_relative_error = [np.abs(self.tm_band_freqs[i+1] - self.tm_band_freqs[i]) for i in range(len(self.tm_band_freqs)-1)]
		tm_absolute_error = [np.abs(self.tm_band_freqs[i] - self.tm_band_freqs[-1]) for i in range(len(self.tm_band_freqs)-1)]
		fig = plt.figure()
		plt.subplot(121)
		plt.title('k='+str(self.k[k_ind]))
		leg = []
		for band in self.bands:
			plt.plot(self.parameters[self.sweep_parameter][self.start:self.stop], np.array(te_relative_error)[self.start:self.stop, k_ind, band - 1])
			plt.plot(self.parameters[self.sweep_parameter][self.start:self.stop], np.array(tm_relative_error)[self.start:self.stop, k_ind, band - 1])
			leg.append('te band ' + str(band))
			leg.append('tm band ' + str(band))
		plt.xlabel(self.sweep_parameter)
		plt.ylabel('Relative Error')
		if log:
			plt.yscale('log') 
		if legend:
			plt.legend(leg, loc='upper right')
		ax = plt.subplot(122)
		plt.title('k='+str(self.k[k_ind]))
		ax.yaxis.tick_right()
		ax.yaxis.set_label_position("right")
		leg = []
		for band in self.bands:
			plt.plot(self.parameters[self.sweep_parameter][self.start:self.stop], np.array(te_absolute_error)[self.start:self.stop, k_ind, band - 1])
			plt.plot(self.parameters[self.sweep_parameter][self.start:self.stop], np.array(tm_absolute_error)[self.start:self.stop, k_ind, band - 1])
			leg.append('te band ' + str(band))
			leg.append('tm band ' + str(band))
		plt.xlabel(self.sweep_parameter)
		plt.ylabel('Absolute Error')
		if log:
			plt.yscale('log')
		if legend:
			plt.legend(leg, loc='upper right')

