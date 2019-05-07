import sys, re
import PyQt5
from PyQt5.QtWidgets import *
from PyQt5 import uic
from PyQt5.QtCore import pyqtSlot, QDate, Qt
from PyQt5.QtGui import QIcon, QPixmap, QFont, QImage

import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.animation import FuncAnimation
import matplotlib.colors
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.patches as mpatches
from scipy import stats

import time
import os
from shutil import copyfile

import MDAnalysis as md
from MDAnalysis.analysis.rdf import InterRDF
import scipy as sp
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster

current_time = lambda : time.time()

dir_principal = os.getcwd()

data_folder = dir_principal + '/Data' #.ui files and temporal .png files
array_folder = os.getcwd() + '/Analysis files'

linear_chains_folder = array_folder + '/Linear chains' #Folder inside analysis files
cluster_folder = array_folder + '/Clusters' #Folder inside analysis file
general_aggregates_folder = array_folder + '/General aggregates' #
isolated_particles_folder = array_folder + '/Isolated particles'
lateral_chains_folder = array_folder + '/Lateral chains'
isolated_aggregates_folder = array_folder + '/Isolated B aggregates'
rdf_folder = array_folder + '/RDF'

if not os.path.exists(data_folder): os.mkdir(data_folder)
if not os.path.exists(array_folder): os.mkdir(array_folder)

if not os.path.exists(linear_chains_folder): os.mkdir(linear_chains_folder)
if not os.path.exists(cluster_folder): os.mkdir(cluster_folder)
if not os.path.exists(lateral_chains_folder): os.mkdir(lateral_chains_folder)
if not os.path.exists(isolated_aggregates_folder): os.mkdir(isolated_aggregates_folder)
if not os.path.exists(isolated_particles_folder): os.mkdir(isolated_particles_folder)
if not os.path.exists(general_aggregates_folder): os.mkdir(general_aggregates_folder)
if not os.path.exists(rdf_folder): os.mkdir(rdf_folder)

os.chdir(data_folder)

class Window(QMainWindow): 
	def __init__(self):
		QMainWindow.__init__(self)
		os.chdir(data_folder)
		uic.loadUi('ChainAnalysis_reduced.ui', self)
		os.chdir(dir_principal)

		self.showMaximized()

		self.manage_files = ManageFiles()

		#Load files
		self.load_topology.clicked.connect(self.obrir_load_files_topology)
		self.load_trajectory.clicked.connect(self.obrir_load_files_trajectory)
		self.load_cluster.clicked.connect(self.obrir_load_files_cluster)
		self.load_linear_chains.clicked.connect(self.obrir_load_files_linear_chains)
		self.load_general_aggregates.clicked.connect(self.obrir_load_files_general_aggregates)
		self.load_isolated_particles.clicked.connect(self.obrir_load_files_isolated_particles)
		self.load_isolated_aggregates.clicked.connect(self.obrir_load_files_isolated_aggregates)
		self.load_rdf.clicked.connect(self.obrir_load_files_rdf)

		self.load_system_param.clicked.connect(self.obrir_load_files_system_param)

		#System parameters from file
		self.read_from_file.clicked.connect(self.obrir_read_from_file)

		#Analysis
		self.linear_chains_finder.clicked.connect(self.obrir_linear_chains_finder)
		self.cluster_finder.clicked.connect(self.obrir_cluster_finder)
		self.general_isolated_finder.clicked.connect(self.obrir_general_isolated_finder)
		self.isolated_aggregates_finder.clicked.connect(self.obrir_isolated_aggregates_finder)
		self.rdf_analysis.clicked.connect(self.obrir_rdf_analysis)

		#Make plots
		self.make_plot.clicked.connect(self.obrir_make_plot)


	#Load files
	def obrir_load_files_topology(self):
		self.manage_files.openFileNameDialog(self.loaded_topology)

	def obrir_load_files_trajectory(self):
		self.manage_files.openFileNameDialog(self.loaded_trajectory)

	def obrir_load_files_cluster(self):
		self.manage_files.openFileNameDialog(self.loaded_cluster)

	def obrir_load_files_linear_chains(self):
		self.manage_files.openFileNameDialog(self.loaded_linear_chains)

	def obrir_load_files_general_aggregates(self):
		self.manage_files.openFileNameDialog(self.loaded_general_aggregates)

	def obrir_load_files_isolated_particles(self):
		self.manage_files.openFileNameDialog(self.loaded_isolated_particles)

	def obrir_load_files_isolated_aggregates(self):
		self.manage_files.openFileNameDialog(self.loaded_isolated_aggregates)

	def obrir_load_files_system_param(self):
		self.manage_files.openFileNameDialog(self.loaded_system_param)

	def obrir_load_files_rdf(self):
		self.manage_files.openFileNameDialog(self.loaded_rdf)


	#System parameters from file
	def obrir_read_from_file(self):
		if self.loaded_system_param.text() != '':
			f = open(self.loaded_system_param.text())

			param = [self.gamma_A, self.gamma_B, self.phi_A, self.phi_B, self.N_A, self.N_B]

			i = 0
			for line in f:
				data = line.split()

				if i == 0 or i == 1:
					param[i].setValue(float(data[-1]))

				elif i == 2 or i == 3:
					param[i].setText(str(data[-1]))

				else:
					param[i].setValue(int(data[-1]))

				i += 1

		else:
			QMessageBox.warning(self, 'Warning!', 'You must load the system parameters file first!')


	#Analysis
	def obrir_linear_chains_finder(self):

		topology = self.loaded_topology.text()
		trajectory = self.loaded_trajectory.text()

		radial_cutoff = self.radial_cutoff.value()
		linear_cutoff = self.linear_cutoff.value()

		init_frame = self.init_frame.value()
		final_frame = self.final_frame.value()
		dt = self.dt.value()

		selection = self.selection.text()

		filename = self.filename_linear_chains.text()

		if filename == '':
			QMessageBox.warning(self, 'Warning!', 'You must name the .npy file to save!')

		elif topology == '':
			QMessageBox.warning(self, 'Warning!', 'You must load the topology file!')

		elif trajectory == '':
			QMessageBox.warning(self, 'Warning!', 'You must load the trajectory file!')

		elif selection == '':
			QMessageBox.warning(self, 'Warning!', 'You must choose a selection!')

		else:
			
			'''
			Find linear aggregates
			'''

			t0 = current_time()

			U = md.Universe(topology, trajectory) #loading trajectory

			selected_particles = U.select_atoms(selection)

			B_in_linear_chains = []

			if final_frame == 'last':
				final_frame = len(U.trajectory) - 1

			total_frames = int(final_frame - init_frame)

			for i in range(init_frame):
				U.trajectory.next()

			for ts in range(0, len(U.trajectory[init_frame:final_frame]), dt): #loop over frames

				N_B_in_linear_chains = 0

				for particle in selected_particles: #Loop over B particles
					index = particle.index

					sel1 = U.select_atoms("cyzone %.2f %.2f %.2f bynum %i" % (radial_cutoff, linear_cutoff, -linear_cutoff, index)) #selection of all particles aggregated with particle of index=index
					neighbours = len(sel1) - 1 #Avoid selfcounting

					if neighbours > 0: #If B particle has some neighbour, it is a B particle in a chain
					    N_B_in_linear_chains += 1

				B_in_linear_chains.append(N_B_in_linear_chains)

				self.progressBar_linear_chains.setValue(int(ts / total_frames * 100))

				#Update trjacetory
				for i in range(dt):
					U.trajectory.next()

			os.chdir(linear_chains_folder)
			np.save(filename, B_in_linear_chains)

			time_used = round((current_time() - t0), 2)

			if time_used < 60:
				self.execution_time.setText(str(time_used) + ' s')	

			else:
				time_used = round(time_used / 60, 2)

				self.execution_time.setText(str(time_used) + ' min')

			QMessageBox.information(self, 'Information', 'Linear Chain analysis finished!')
			self.progressBar_linear_chains.setValue(0)	

	def obrir_cluster_finder(self):

		topology = self.loaded_topology.text()
		trajectory = self.loaded_trajectory.text()

		radial_cutoff = self.radial_cutoff.value()
		linear_cutoff = self.linear_cutoff.value()

		init_frame = self.init_frame.value()
		final_frame = self.final_frame.value()
		dt = self.dt.value()

		selection = self.selection.text()

		filename = self.filename_cluster.text()

		if filename == '':
			QMessageBox.warning(self, 'Warning!', 'You must name the .npy file to save!')

		elif topology == '':
			QMessageBox.warning(self, 'Warning!', 'You must load the topology file!')

		elif trajectory == '':
			QMessageBox.warning(self, 'Warning!', 'You must load the trajectory file!')

		elif selection == '':
			QMessageBox.warning(self, 'Warning!', 'You must choose a selection!')

		else:

			t0 = current_time()

			U = md.Universe(topology, trajectory) #loading trajectory

			sel1 = U.select_atoms(selection) #selection

			aggregates_T = []
			total_frames = int(final_frame - init_frame)

			for i in range(init_frame):
				U.trajectory.next()

			for ts in range(0, len(U.trajectory[init_frame:final_frame]), dt):
				self.progressBar_cluster.setValue(int((ts / total_frames) * 100))

				distan = md.lib.distances.self_distance_array(sel1.positions,box=U.dimensions)#distance matrix
				clust = sp.cluster.hierarchy.linkage(distan,'single') #clustering using complete linkage method

				clustdist = fcluster(clust, linear_cutoff, criterion='distance')#clasifica los clusters con distancias menores a cutoff
				nclust = np.amax(clustdist) #number of clusters 

				#distribucion de clusters
				binss = nclust

				clustsize = np.histogram(clustdist, bins=binss)

				aggregates = []

				#Append the aggregates in each time step
				for item in clustsize[0]:
					if item > 1:
						aggregates.append(item)

				aggregates_T.append(aggregates)

				#Update trajectory
				for i in range(dt):
					U.trajectory.next()

			aggregates_T = np.array(aggregates_T)

			os.chdir(cluster_folder)
			np.save(filename, aggregates_T)

			time_used = round((current_time() - t0), 2)

			if time_used < 60:
				self.execution_time.setText(str(time_used) + ' s')	

			else:
				time_used = round(time_used / 60, 2)

				self.execution_time.setText(str(time_used) + ' min')

			QMessageBox.information(self, 'Information', 'Cluster analysis finished!')
			self.progressBar_cluster.setValue(0)

	def obrir_lateral_chains_finder(self):
		'''Calculate linear chains again with a higher radial cutoff in order to
		   make the difference later and obtain number of particles aggregated laterally 
		'''

		topology = self.loaded_topology.text()
		trajectory = self.loaded_trajectory.text()

		radial_cutoff = self.radial_cutoff.value()
		linear_cutoff = self.linear_cutoff.value()

		init_frame = self.init_frame.value()
		final_frame = self.final_frame.value()
		dt = self.dt.value()

		selection = self.selection.text()

		filename = self.filename_lateral_chains.text()

		if filename == '':
			QMessageBox.warning(self, 'Warning!', 'You must name the .npy file to save!')

		elif topology == '':
			QMessageBox.warning(self, 'Warning!', 'You must load the topology file!')

		elif trajectory == '':
			QMessageBox.warning(self, 'Warning!', 'You must load the trajectory file!')

		elif selection == '':
			QMessageBox.warning(self, 'Warning!', 'You must choose a selection!')

		else:
			
			'''
			Find linear aggregates
			'''

			t0 = current_time()

			U = md.Universe(topology, trajectory) #loading trajectory

			selected_particles = U.select_atoms(selection)

			B_in_linear_chains = []

			if final_frame == 'last':
				final_frame = len(U.trajectory) - 1

			total_frames = int(final_frame - init_frame)

			for i in range(init_frame):
				U.trajectory.next()

			for ts in range(0, len(U.trajectory[init_frame:final_frame]), dt): #loop over frames

				N_B_in_linear_chains = 0

				for particle in selected_particles: #Loop over B particles
					index = particle.index

					sel1 = U.select_atoms("cyzone %.2f %.2f %.2f bynum %i" % (radial_cutoff, linear_cutoff, -linear_cutoff, index)) #selection of all particles aggregated with particle of index=index
					neighbours = len(sel1) - 1 #Avoid selfcounting

					if neighbours > 0: #If B particle has some neighbour, it is a B particle in a chain
					    N_B_in_linear_chains += 1

				B_in_linear_chains.append(N_B_in_linear_chains)

				self.progressBar_linear_chains.setValue(int(ts / total_frames * 100))

				#Update trjacetory
				for i in range(dt):
					U.trajectory.next()

			os.chdir(linear_chains_folder)
			np.save(filename, B_in_linear_chains)

			time_used = round((current_time() - t0), 2)

			if time_used < 60:
				self.execution_time.setText(str(time_used) + ' s')	

			else:
				time_used = round(time_used / 60, 2)

				self.execution_time.setText(str(time_used) + ' min')

			QMessageBox.information(self, 'Information', 'Lateral Chain analysis finished!')
			self.progressBar_lateral_chains.setValue(0)	

	def obrir_general_isolated_finder(self):
		topology = self.loaded_topology.text()
		trajectory = self.loaded_trajectory.text()

		radial_cutoff = self.radial_cutoff.value()

		init_frame = self.init_frame.value()
		final_frame = self.final_frame.value()
		dt = self.dt.value()

		selection = self.selection.text()

		filenames = self.filename_general_isolated_finder.text().split()

		if len(filenames) < 2:
			QMessageBox.warning(self, 'Warning!', 'You must choose 2 filenames to save general and isolated particles (in order)')

		elif topology == '':
			QMessageBox.warning(self, 'Warning!', 'You must load the topology file!')

		elif trajectory == '':
			QMessageBox.warning(self, 'Warning!', 'You must load the trajectory file!')

		elif selection == '':
			QMessageBox.warning(self, 'Warning!', 'You must choose a selection!')

		else:

			filename = filenames[0]
			filename_2 = filenames[1]
			
			'''Find general aggregates and isolated particles of the current selection'''

			t0 = current_time()

			U = md.Universe(topology, trajectory) #loading trajectory

			selected_particles = U.select_atoms(selection) #User's selection

			isolated_particles = [] #Number of isolated particles in time
			aggregated_particles = [] #Number of aggregated particles in time

			for i in range(init_frame):
				U.trajectory.next()

			total_frames = int(final_frame - init_frame)

			for ts in range(0, len(U.trajectory[init_frame:final_frame]), dt): #loop over frames

				isolated = 0
				aggregated = 0

				for particle in selected_particles: #Loop over B particles
					index = particle.index

					sel1 = U.select_atoms("sphzone %.2f bynum %i" % (radial_cutoff, index)) #selection of all particles aggregated with particle of index=index
					neighbours = len(sel1) - 1 #Avoid selfcounting

					if neighbours == 0: #If B particle has some neighbour, it is a B particle in a chain
						isolated +=1

					else:
						aggregated += 1

				isolated_particles.append(isolated)
				aggregated_particles.append(aggregated)

				#Update trjacetory
				for i in range(dt):
					U.trajectory.next()

				self.progressBar_general_isolated.setValue(int(ts / total_frames * 100))

			os.chdir(isolated_particles_folder)
			np.save(filename, isolated_particles)

			os.chdir(general_aggregates_folder)
			np.save(filename_2, aggregated_particles)

			time_used = round((current_time() - t0), 2)

			if time_used < 60:
				self.execution_time.setText(str(time_used) + ' s')	

			else:
				time_used = round(time_used / 60, 2)

				self.execution_time.setText(str(time_used) + ' min')

			QMessageBox.information(self, 'Information', 'Analysis finished!')
			self.progressBar_linear_chains.setValue(0)

	def obrir_isolated_aggregates_finder(self):

		topology = self.loaded_topology.text()
		trajectory = self.loaded_trajectory.text()

		linear_cutoff = self.linear_cutoff.value()

		init_frame = self.init_frame.value()
		final_frame = self.final_frame.value()
		dt = self.dt.value()

		selection = self.selection.text()

		filename = self.filename_isolated_aggregates.text()

		total_frames = int(final_frame - init_frame)

		if filename == '':
			QMessageBox.warning(self, 'Warning!', 'You must name the .npy file to save!')

		elif topology == '':
			QMessageBox.warning(self, 'Warning!', 'You must load the topology file!')

		elif trajectory == '':
			QMessageBox.warning(self, 'Warning!', 'You must load the trajectory file!')

		elif selection == '':
			QMessageBox.warning(self, 'Warning!', 'You must choose a selection!')

		else:
			t0 = current_time()

			U = md.Universe(topology, trajectory) #loading trajectory

			selected_particles = U.select_atoms(selection)

			isolated_B_agg = []

			for i in range(init_frame):
				U.trajectory.next()

			for ts in range(0, len(U.trajectory[init_frame:final_frame]), dt): #loop over frames

				self.progressBar_isolated_aggregates.setValue(int(ts / total_frames * 100))

				B_agg = []

				for particle in selected_particles: #Loop over B particles
					index = int(particle.index) + 1

					sel1 = U.select_atoms("sphlayer 0.5 %.2f bynum %i" % (linear_cutoff, index)) #selection of all particles aggregated with particle of index=index

					if len(sel1) != 0:

						#Check if it's aggregated with some A
						A_agg = False
							
						for part in sel1:
							if int(part.type) == 1: #The current B particle is aggregated with another A particle
								A_agg = True
							
						if A_agg == False: #If it's only aggregated with B, append to the list

							B_agg.append(index)

				lenght_not_stable = True

				while lenght_not_stable:

					l1 = len(B_agg)

					for index in B_agg:

						sel = U.select_atoms("sphlayer 0.5 %.2f bynum %i" % (linear_cutoff, index)) #Particles aggregated with the current B part

						for particle in sel:
							current_index = int(particle.index) + 1 #index of one of the aggregated particles

							if not current_index in B_agg: #If some of the aggregated particles is not in the white list it's not a isolated_B_agg
								try:
									B_agg.remove(index)
								except:
									pass

					l2 = len(B_agg)

					L = l2 - l1

					if L == 0:
						lenght_not_stable = False
				

				isolated_B_agg.append(len(B_agg)) #B with B's isolated (not in chains or aggregates thah contain A)

				for i in range(dt):
					U.trajectory.next()


			os.chdir(isolated_aggregates_folder)
			np.save(filename, isolated_B_agg)

			time_used = round((current_time() - t0), 2)

			if time_used < 60:
				self.execution_time.setText(str(time_used) + ' s')	

			else:
				time_used = round(time_used / 60, 2)

				self.execution_time.setText(str(time_used) + ' min')

			QMessageBox.information(self, 'Information', 'Analysis finished!')
			self.progressBar_isolated_aggregates.setValue(0)

	def obrir_rdf_analysis(self):

		topology = self.loaded_topology.text()
		trajectory = self.loaded_trajectory.text()

		filename = self.filename_rdf.text()

		init_frame = self.init_frame.value()
		final_frame = self.final_frame.value()
		dt = self.dt.value()

		r_max = self.r_max.value()
		dr = self.dr.value()

		bins = np.arange(0, r_max, dr)

		selection = self.selection.text()
		selection_2 = self.selection_2.text()

		if topology == '':
			QMessageBox.warning(self, 'Warning!', 'You must load the topology file!')

		elif trajectory == '':
			QMessageBox.warning(self, 'Warning!', 'You must load the trajectory file!')

		elif filename == '':
			QMessageBox.warning(self, 'Warning!', 'You must choose a filename for the analysed data output!')

		elif selection == '' or selection_2 == '':
			QMessageBox.warning(self, 'Warning!', 'You must choose a two selections') 

		else:

			t0 = current_time()

			self.progress_rdf.setText('Calculating...')

			U = md.Universe(topology, trajectory) #loading trajectory

			B_particles = U.select_atoms(selection_2)
			A_particles = U.select_atoms(selection)

			rdf = InterRDF(A_particles, B_particles, nbins=bins, range=(0.0, r_max), exclusion_block=(1,1), verbose=False)
			rdf.run(start=init_frame, stop=final_frame, step=dt)

			self.progress_rdf.setText('Done!')

			os.chdir(rdf_folder)
			np.save(filename, [rdf.bins, rdf.rdf])

			time_used = round((current_time() - t0), 2)

			QMessageBox.information(self, 'Information', 'RDF analysis finished!')

			if time_used < 60:
					self.execution_time.setText(str(time_used) + ' s')	

			else:
				time_used = round(time_used / 60, 2)

				self.execution_time.setText(str(time_used) + ' min')


	#Choose plot to make
	def obrir_make_plot(self):
		if self.plot_type.currentText() == 'Particles aggregated linearly':
			self.make_plot_B_linear_chains()

		elif self.plot_type.currentText() == 'Particles aggregated laterally':
			self.make_plot_B_lateral_chains()

		elif self.plot_type.currentText() == 'Number of clusters':
			self.make_plot_N_clusters()

		elif self.plot_type.currentText() == 'Avg cluster length':
			self.make_plot_avg_clusters_length()

		elif self.plot_type.currentText() == 'Isolated particles':
			self.make_plot_isolated_particles()

		elif self.plot_type.currentText() == 'Aggregated particles':
			self.make_plot_general_aggregates()

		elif self.plot_type.currentText() == 'Isolated aggregates':
			self.make_plot_isolated_aggregates()

		elif self.plot_type.currentText() == 'q linear aggregates':
			self.make_plot_q_linear_chains()

		elif self.plot_type.currentText() == 'q general aggregates':
			self.make_plot_q_lateral_chains()

		elif self.plot_type.currentText() == 'Radial Distribution Function':
			self.make_plot_rdf()


	#Make plots
	def make_plot_N_clusters(self):
		#System parameters
		gamma_A = self.gamma_A.value()
		gamma_B = self.gamma_B.value()

		phi_A = float(self.phi_A.text())
		phi_B = float(self.phi_B.text())

		N_A = self.N_A.value()
		N_B = self.N_B.value()

		#Analysis parameters
		init_frame = self.init_frame_plot.value()
		final_frame = self.final_frame_plot.value()
		dt = self.dt_plot.value()
		eq_time = self.eq_time.value()

		title = self.title_plot.text()
		xlabel = self.xlabel_plot.text()
		ylabel = self.ylabel_plot.text()
		labels = self.labels_plot.text()

		linestyles = self.ls_plot.text()
		colors = self.colors_plot.text()

		plot_avg = self.average_plot.isChecked() #Checkbox State (true or false)
		last_frame = self.last_frame.isChecked()

		array_filename = self.loaded_cluster.text()

		default = False

		if labels == '' or colors == '' or linestyles == '':
			default = True

		elif self.default_plot_settings.isChecked():
			default = True

		if array_filename != '':

			plot = MakePlot(init_frame, final_frame, dt, eq_time, title, xlabel, ylabel, labels, linestyles, colors, plot_avg, last_frame, array_filename, gamma_A, gamma_B, phi_A, phi_B, N_A, N_B)

			t, N_aggregates_t, avg_N_aggregates, avg_N, Average_lenght_aggregates_t, avg_avg_length, avg_l, last_frame = plot.aggregate_analysis() #data analysis to make the plot

			if final_frame > last_frame:
				QMessageBox.warning(self, 'Warning!', 'Final frame max value exceeded. Already set to max value.')

			self.final_frame_plot.setValue(last_frame)

			plot.N_clusters(t, N_aggregates_t, avg_N_aggregates, avg_N, 'N_clusters.png', default) #Saves the plot as .png file

		else:
			QMessageBox.warning(self, 'Warning!', 'You must load the cluster data file!')

	def make_plot_avg_clusters_length(self):
		#System parameters
		gamma_A = self.gamma_A.value()
		gamma_B = self.gamma_B.value()

		phi_A = float(self.phi_A.text())
		phi_B = float(self.phi_B.text())

		N_A = self.N_A.value()
		N_B = self.N_B.value()

		#Analysis parameters
		init_frame = self.init_frame_plot.value()
		final_frame = self.final_frame_plot.value()
		dt = self.dt_plot.value()
		eq_time = self.eq_time.value()

		title = self.title_plot.text()
		xlabel = self.xlabel_plot.text()
		ylabel = self.ylabel_plot.text()
		labels = self.labels_plot.text()

		linestyles = self.ls_plot.text()
		colors = self.colors_plot.text()

		plot_avg = self.average_plot.isChecked() #Checkbox State (true or false)
		last_frame = self.last_frame.isChecked()

		array_filename = self.loaded_cluster.text()

		default = False

		if labels == '' or colors == '' or linestyles == '':
			default = True

		elif self.default_plot_settings.isChecked():
			default = True

		if array_filename != '':

			plot = MakePlot(init_frame, final_frame, dt, eq_time, title, xlabel, ylabel, labels, linestyles, colors, plot_avg, last_frame, array_filename, gamma_A, gamma_B, phi_A, phi_B, N_A, N_B)

			t, N_aggregates_t, avg_N_aggregates, avg_N, Average_lenght_aggregates_t, avg_avg_length, avg_l, last_frame = plot.aggregate_analysis() #data analysis to make the plot

			if final_frame > last_frame:
				QMessageBox.warning(self, 'Warning!', 'Final frame max value exceeded. Already set to max value.')

			self.final_frame_plot.setValue(last_frame)

			plot.avg_cluster_length(t, Average_lenght_aggregates_t, avg_avg_length, avg_l, 'avg_length_clusters.png', default) #Saves the plot as .png file

		else:
			QMessageBox.warning(self, 'Warning!', 'You must load the cluster data file!')

	def make_plot_B_linear_chains(self):
		#System parameters
		gamma_A = self.gamma_A.value()
		gamma_B = self.gamma_B.value()

		phi_A = float(self.phi_A.text())
		phi_B = float(self.phi_B.text())

		N_A = self.N_A.value()
		N_B = self.N_B.value()

		#Analysis parameters
		init_frame = self.init_frame_plot.value()
		final_frame = self.final_frame_plot.value()
		dt = self.dt_plot.value()
		eq_time = self.eq_time.value()

		title = self.title_plot.text()
		xlabel = self.xlabel_plot.text()
		ylabel = self.ylabel_plot.text()
		labels = self.labels_plot.text()

		linestyles = self.ls_plot.text()
		colors = self.colors_plot.text()

		plot_avg = self.average_plot.isChecked() #Checkbox State (true or false)
		last_frame = self.last_frame.isChecked()

		array_filename = self.loaded_linear_chains.text()

		default = False

		if labels == '' or colors == '' or linestyles == '':
			default = True

		elif self.default_plot_settings.isChecked():
			default = True

		if array_filename != '':

			plot = MakePlot(init_frame, final_frame, dt, eq_time, title, xlabel, ylabel, labels, linestyles, colors, plot_avg, last_frame, array_filename, gamma_A, gamma_B, phi_A, phi_B, N_A, N_B)

			t, B_in_chains, avg_B_in_chains, qs, avg_q, last_frame = plot.linear_chains_statistics() #data analysis to make the plot

			if final_frame > last_frame:
				QMessageBox.warning(self, 'Warning!', 'Final frame max value exceeded. Already set to max value.')

			self.final_frame_plot.setValue(last_frame)

			plot.N_B_in_linear_chains(t, B_in_chains, avg_B_in_chains, 'B_linear_chains.png', default) #Saves the plot as .png file

		else:
			QMessageBox.warning(self, 'Warning!', 'You must load the linear chains data file!')

	def make_plot_q_linear_chains(self):
		#System parameters
		gamma_A = self.gamma_A.value()
		gamma_B = self.gamma_B.value()

		phi_A = float(self.phi_A.text())
		phi_B = float(self.phi_B.text())

		N_A = self.N_A.value()
		N_B = self.N_B.value()

		#Analysis parameters
		init_frame = self.init_frame_plot.value()
		final_frame = self.final_frame_plot.value()
		dt = self.dt_plot.value()
		eq_time = self.eq_time.value()

		title = self.title_plot.text()
		xlabel = self.xlabel_plot.text()
		ylabel = self.ylabel_plot.text()
		labels = self.labels_plot.text()

		linestyles = self.ls_plot.text()
		colors = self.colors_plot.text()

		plot_avg = self.average_plot.isChecked() #Checkbox State (true or false)
		last_frame = self.last_frame.isChecked()

		array_filename = self.loaded_linear_chains.text()

		default = False

		if labels == '' or colors == '' or linestyles == '':
			default = True

		elif self.default_plot_settings.isChecked():
			default = True

		if array_filename != '':

			plot = MakePlot(init_frame, final_frame, dt, eq_time, title, xlabel, ylabel, labels, linestyles, colors, plot_avg, last_frame, array_filename, gamma_A, gamma_B, phi_A, phi_B, N_A, N_B)

			t, B_in_chains, avg_B_in_chains, qs, avg_q, last_frame = plot.linear_chains_statistics() #data analysis to make the plot

			if final_frame > last_frame:
				QMessageBox.warning(self, 'Warning!', 'Final frame max value exceeded. Already set to max value.')

			self.final_frame_plot.setValue(last_frame)

			plot.q_linear_chains(t, qs, avg_q, 'q_linear_chains.png', default) #Saves the plot as .png file

		else:
			QMessageBox.warning(self, 'Warning!', 'You must load the linear chains data file!')

	def make_plot_B_lateral_chains(self):
		#System parameters
		gamma_A = self.gamma_A.value()
		gamma_B = self.gamma_B.value()

		phi_A = float(self.phi_A.text())
		phi_B = float(self.phi_B.text())

		N_A = self.N_A.value()
		N_B = self.N_B.value()

		#Analysis parameters
		init_frame = self.init_frame_plot.value()
		final_frame = self.final_frame_plot.value()
		dt = self.dt_plot.value()
		eq_time = self.eq_time.value()

		title = self.title_plot.text()
		xlabel = self.xlabel_plot.text()
		ylabel = self.ylabel_plot.text()
		labels = self.labels_plot.text()

		linestyles = self.ls_plot.text()
		colors = self.colors_plot.text()

		plot_avg = self.average_plot.isChecked() #Checkbox State (true or false)
		last_frame = self.last_frame.isChecked()

		filename1 = self.loaded_linear_chains.text()
		array_filename = self.loaded_general_aggregates.text()

		default = False

		if labels == '' or colors == '' or linestyles == '':
			default = True

		elif self.default_plot_settings.isChecked():
			default = True

		if array_filename != '' and filename1 != '':

			plot = MakePlot(init_frame, final_frame, dt, eq_time, title, xlabel, ylabel, labels, linestyles, colors, plot_avg, last_frame, array_filename, gamma_A, gamma_B, phi_A, phi_B, N_A, N_B)

			t, B_in_lateral_chains, avg_B_in_lateral_chains, qs, avg_q, last_frame, different_length = plot.lateral_chains_statistics(filename1, array_filename) #data analysis to make the plot

			if different_length == True:
				QMessageBox.warning(self, 'Warning!', 'The arrays loaded don\'t have the same length! Using the shortest length for the analysis.')

			if final_frame > last_frame:
				QMessageBox.warning(self, 'Warning!', 'Final frame max value exceeded. Already set to max value.')

			self.final_frame_plot.setValue(last_frame)

			plot.N_B_in_linear_chains(t, B_in_lateral_chains, avg_B_in_lateral_chains, 'B_lateral_chains.png', default) #Saves the plot as .png file

		else:
			QMessageBox.warning(self, 'Warning!', 'You must load the linear chains and general aggregates data files!')

	def make_plot_q_lateral_chains(self):
		#System parameters
		gamma_A = self.gamma_A.value()
		gamma_B = self.gamma_B.value()

		phi_A = float(self.phi_A.text())
		phi_B = float(self.phi_B.text())

		N_A = self.N_A.value()
		N_B = self.N_B.value()

		#Analysis parameters
		init_frame = self.init_frame_plot.value()
		final_frame = self.final_frame_plot.value()
		dt = self.dt_plot.value()
		eq_time = self.eq_time.value()

		title = self.title_plot.text()
		xlabel = self.xlabel_plot.text()
		ylabel = self.ylabel_plot.text()
		labels = self.labels_plot.text()

		linestyles = self.ls_plot.text()
		colors = self.colors_plot.text()

		plot_avg = self.average_plot.isChecked() #Checkbox State (true or false)
		last_frame = self.last_frame.isChecked()

		filename1 = self.loaded_linear_chains.text()
		array_filename = self.loaded_general_aggregates.text()

		default = False

		if labels == '' or colors == '' or linestyles == '':
			default = True

		elif self.default_plot_settings.isChecked():
			default = True

		if array_filename != '' and filename1 != '':

			plot = MakePlot(init_frame, final_frame, dt, eq_time, title, xlabel, ylabel, labels, linestyles, colors, plot_avg, last_frame, array_filename, gamma_A, gamma_B, phi_A, phi_B, N_A, N_B)

			t, B_in_lateral_chains, avg_B_in_lateral_chains, qs, avg_q, last_frame, different_length = plot.lateral_chains_statistics(filename1, array_filename) #data analysis to make the plot

			if different_length == True:
				QMessageBox.warning(self, 'Warning!', 'The arrays loaded don\'t have the same length! Using the shortest length for the analysis.')

			if final_frame > last_frame:
				QMessageBox.warning(self, 'Warning!', 'Final frame max value exceeded. Already set to max value.')

			self.final_frame_plot.setValue(last_frame)

			plot.q_linear_chains(t, qs, avg_q, 'q_linear_lateral_chains.png', default) 

		else:
			QMessageBox.warning(self, 'Warning!', 'You must load the linear chains and general aggregates data files!')

	def make_plot_general_aggregates(self):
		#System parameters
		gamma_A = self.gamma_A.value()
		gamma_B = self.gamma_B.value()

		phi_A = float(self.phi_A.text())
		phi_B = float(self.phi_B.text())

		N_A = self.N_A.value()
		N_B = self.N_B.value()

		#Analysis parameters
		init_frame = self.init_frame_plot.value()
		final_frame = self.final_frame_plot.value()
		dt = self.dt_plot.value()
		eq_time = self.eq_time.value()

		title = self.title_plot.text()
		xlabel = self.xlabel_plot.text()
		ylabel = self.ylabel_plot.text()
		labels = self.labels_plot.text()

		linestyles = self.ls_plot.text()
		colors = self.colors_plot.text()

		plot_avg = self.average_plot.isChecked() #Checkbox State (true or false)
		last_frame = self.last_frame.isChecked()

		array_filename = self.loaded_general_aggregates.text()

		default = False

		if labels == '' or colors == '' or linestyles == '':
			default = True

		elif self.default_plot_settings.isChecked():
			default = True

		if array_filename != '':

			plot = MakePlot(init_frame, final_frame, dt, eq_time, title, xlabel, ylabel, labels, linestyles, colors, plot_avg, last_frame, array_filename, gamma_A, gamma_B, phi_A, phi_B, N_A, N_B)

			t, N_aggregates, avg, last_frame = plot.general_aggregates_analysis() #data analysis to make the plot

			if final_frame > last_frame:
				QMessageBox.warning(self, 'Warning!', 'Final frame max value exceeded. Already set to max value.')

			self.final_frame_plot.setValue(last_frame)

			plot.N_general_aggregates(t, N_aggregates, avg, 'general_aggregates.png', default) #Saves the plot as .png file

		else:
			QMessageBox.warning(self, 'Warning!', 'You must load the general aggregates data file!')

	def make_plot_isolated_particles(self):
		#System parameters
		gamma_A = self.gamma_A.value()
		gamma_B = self.gamma_B.value()

		phi_A = float(self.phi_A.text())
		phi_B = float(self.phi_B.text())

		N_A = self.N_A.value()
		N_B = self.N_B.value()

		#Analysis parameters
		init_frame = self.init_frame_plot.value()
		final_frame = self.final_frame_plot.value()
		dt = self.dt_plot.value()
		eq_time = self.eq_time.value()

		title = self.title_plot.text()
		xlabel = self.xlabel_plot.text()
		ylabel = self.ylabel_plot.text()
		labels = self.labels_plot.text()

		linestyles = self.ls_plot.text()
		colors = self.colors_plot.text()

		plot_avg = self.average_plot.isChecked() #Checkbox State (true or false)
		last_frame = self.last_frame.isChecked()

		array_filename = self.loaded_isolated_particles.text()

		default = False

		if labels == '' or colors == '' or linestyles == '':
			default = True

		elif self.default_plot_settings.isChecked():
			default = True

		if array_filename != '':

			plot = MakePlot(init_frame, final_frame, dt, eq_time, title, xlabel, ylabel, labels, linestyles, colors, plot_avg, last_frame, array_filename, gamma_A, gamma_B, phi_A, phi_B, N_A, N_B)

			t, isolated_particles, avg, last_frame = plot.isolated_particles_analysis() #data analysis to make the plot

			if final_frame > last_frame:
				QMessageBox.warning(self, 'Warning!', 'Final frame max value exceeded. Already set to max value.')

			self.final_frame_plot.setValue(last_frame)

			plot.N_isolated_particles(t, isolated_particles, avg, 'isolated_particles.png', default) #Saves the plot as .png file

		else:
			QMessageBox.warning(self, 'Warning!', 'You must load the isolated particles data file!')

	def make_plot_isolated_aggregates(self):
		#System parameters
		gamma_A = self.gamma_A.value()
		gamma_B = self.gamma_B.value()

		phi_A = float(self.phi_A.text())
		phi_B = float(self.phi_B.text())

		N_A = self.N_A.value()
		N_B = self.N_B.value()

		#Analysis parameters
		init_frame = self.init_frame_plot.value()
		final_frame = self.final_frame_plot.value()
		dt = self.dt_plot.value()
		eq_time = self.eq_time.value()

		title = self.title_plot.text()
		xlabel = self.xlabel_plot.text()
		ylabel = self.ylabel_plot.text()
		labels = self.labels_plot.text()

		linestyles = self.ls_plot.text()
		colors = self.colors_plot.text()

		plot_avg = self.average_plot.isChecked() #Checkbox State (true or false)
		last_frame = self.last_frame.isChecked()

		array_filename = self.loaded_isolated_aggregates.text()

		default = False

		if labels == '' or colors == '' or linestyles == '':
			default = True

		elif self.default_plot_settings.isChecked():
			default = True

		if array_filename != '':

			plot = MakePlot(init_frame, final_frame, dt, eq_time, title, xlabel, ylabel, labels, linestyles, colors, plot_avg, last_frame, array_filename, gamma_A, gamma_B, phi_A, phi_B, N_A, N_B)

			t, isolated_aggregates, avg, last_frame = plot.isolated_aggregates_analysis() #data analysis to make the plot

			if final_frame > last_frame:
				QMessageBox.warning(self, 'Warning!', 'Final frame max value exceeded. Already set to max value.')

			self.final_frame_plot.setValue(last_frame)

			plot.N_isolated_particles(t, isolated_aggregates, avg, 'isolated_particles.png', default)

		else:
			QMessageBox.warning(self, 'Warning!', 'You must load the isolated particles data file!')

	def make_plot_rdf(self):
		#System parameters
		gamma_A = self.gamma_A.value()
		gamma_B = self.gamma_B.value()

		phi_A = float(self.phi_A.text())
		phi_B = float(self.phi_B.text())

		N_A = self.N_A.value()
		N_B = self.N_B.value()

		#Analysis parameters
		init_frame = self.init_frame_plot.value()
		final_frame = self.final_frame_plot.value()
		dt = self.dt_plot.value()
		eq_time = self.eq_time.value()

		title = self.title_plot.text()
		xlabel = self.xlabel_plot.text()
		ylabel = self.ylabel_plot.text()
		labels = self.labels_plot.text()

		linestyles = self.ls_plot.text()
		colors = self.colors_plot.text()

		plot_avg = self.average_plot.isChecked() #Checkbox State (true or false)
		last_frame = self.last_frame.isChecked()

		array_filename = self.loaded_rdf.text()

		default = False

		if labels == '' or colors == '' or linestyles == '':
			default = True

		elif self.default_plot_settings.isChecked():
			default = True

		if array_filename != '':

			plot = MakePlot(init_frame, final_frame, dt, eq_time, title, xlabel, ylabel, labels, linestyles, colors, plot_avg, last_frame, array_filename, gamma_A, gamma_B, phi_A, phi_B, N_A, N_B)

			plot.plot_rdf(array_filename, 'rdf.png', default) #Saves the plot as .png file

		else:
			QMessageBox.warning(self, 'Warning!', 'You must load the rdf data file!')


	#close event
	def closeEvent(self, event):
		os.chdir(data_folder)
		result = QMessageBox.question(self, 'Leaving...','Do you want to exit?', QMessageBox.Yes | QMessageBox.No)
		if result == QMessageBox.Yes:
			event.accept()	
		else:
			event.ignore()

class ManageFiles(QFileDialog):
	def __init__(self):
		QFileDialog.__init__(self)

		self.title = 'Save files'
		self.left = 10
		self.top = 10
		self.width = 640
		self.height = 400 

		self.initUI()

	def initUI(self):
		self.setWindowTitle(self.title)
		self.setGeometry(self.left, self.top, self.width, self.height)

	def saveFileDialog(self, name):
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog

		fileName, _ = QFileDialog.getSaveFileName(self, 'Save files') 

		if fileName:
			os.chdir(data_folder)
			if os.path.exists('%s.png' % name): copyfile('%s.png' % name, fileName + '.png')
			else: QMessageBox.warning(self, 'Warning!', 'The plot doesn\'t exist!') 

	def openFileNameDialog(self, name):
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog

		fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)

		if fileName:
			name.setText(fileName)

'''
	Chain analysis object class. Can not be imported in the graphical interface class because of ProgressBar. Left here for possible documentation
'''

class ChainAnalysis(object):
	'''
	Topology: .pdb file with the initial coordinates state 

	Trajectory: .dcd file with all the coordinates for all the frames
	'''

	def __init__(self, topology, trajectory, radial_cutoff, linear_cutoff, init_frame, final_frame, dt, selection):
		self.topology = topology
		self.trajectory = trajectory

		self.radial_cutoff = radial_cutoff
		self.linear_cutoff = linear_cutoff

		self.init_frame = init_frame
		self.final_frame = final_frame
		self.dt = dt

		self.selection = selection

	def linear_chains_finder(self, filename):
		'''
		Find linear aggregates
		'''

		t0 = current_time()

		U = md.Universe(self.topology, self.trajectory) #loading trajectory

		selected_particles = self.selection

		B_in_linear_chains = []

		if self.final_frame == 'last':
			self.final_frame = len(U.trajectory) - 1

		total_frames = self.final_frame - self.init_frame

		for ts in range(len(U.trajectory[self.init_frame: self.final_frame])): #loop over frames

			N_B_in_linear_chains = 0

			for particle in selected_particles: #Loop over B particles
				index = particle.index

				sel1 = U.select_atoms("cyzone %.2f %.2f %.2f bynum %i" % (self.radial_cutoff, self.linear_cutoff, -self.linear_cutoff, index)) #selection of all particles aggregated with particle of index=index
				neighbours = len(sel1) - 1 #Avoid selfcounting

				if neighbours > 0: #If B particle has some neighbour, it is a B particle in a chain
				    N_B_in_linear_chains += 1

			B_in_linear_chains.append(N_B_in_linear_chains)

			#Update trajectory
			for i in range(dt):
				U.trajectory.next()

		os.chdir(array_folder)
		np.save(filename, B_in_linear_chains)

		time_used = (current_time() - t0) #About 25 min for 

	def clusters_in_time(self, filename): 

		t0 = current_time()

		U = md.Universe(self.topology, self.trajectory) #loading trajectory

		sel1 = U.select_atoms(self.selection) #selection

		aggregates_T = []

		for ts in range(len(U.trajectory[self.init_frame: self.final_frame])):
			distan = md.lib.distances.self_distance_array(sel1.positions,box=U.dimensions)#distance matrix
			clust = sp.cluster.hierarchy.linkage(distan,'single') #clustering using complete linkage method

			clustdist = fcluster(clust, self.linear_cutoff, criterion='distance')#clasifica los clusters con distancias menores a cutoff
			nclust = np.amax(clustdist) #number of clusters (also clustdists[-1] ???)

			#distribucion de clusters
			binss = nclust

			clustsize = np.histogram(clustdist, bins=binss)

			aggregates = []

			#Append the aggregates in each time step
			for item in clustsize[0]:
				if item > 1:
					aggregates.append(item)

			aggregates_T.append(aggregates)

			#Update trajectory
			for i in range(dt):
				U.trajectory.next()

		delta_T = current_time() - t0

		aggregates_T = np.array(aggregates_T)
		np.save(filename, aggregates_T) #About 5 min for 900 frames

	def isolated_particles_finder(self, filename, filename_2):
		'''
		Find isolated particles
		'''

		t0 = current_time()

		U = md.Universe(self.topology, self.trajectory) #loading trajectory

		selected_particles = U.select_atoms(self.selection)

		isolated_particles = [] #Number of isolated particles in time
		aggregated_particles = [] #Number of aggregated particles in time

		if self.final_frame == 'last':
			self.final_frame = len(U.trajectory) - 1

		total_frames = self.final_frame - self.init_frame

		for ts in range(0, len(U.trajectory[self.init_frame: self.final_frame]), self.dt): #loop over frames
			print(ts)
			isolated = 0
			aggregated = 0

			for particle in selected_particles: #Loop over selected particles
				index = particle.index

				sel1 = U.select_atoms("sphzone %.2f bynum %i" % (self.radial_cutoff, index)) #selection of all particles aggregated with particle of index=index
				neighbours = len(sel1) - 1 #Avoid selfcounting

				if neighbours == 0: #If B particle has some neighbour, it is a B particle in a chain
				    isolated +=1

				else:
					aggregated += 1

			isolated_particles.append(isolated)
			aggregated_particles.append(aggregated)

			#Update trajectory
			for i in range(self.dt):
				U.trajectory.next()

		os.chdir(array_folder)
		np.save(filename, isolated_particles)
		np.save(filename_2, aggregated_particles)

		time_used = (current_time() - t0)

class MakePlot(object):
	def __init__(self, init_frame, final_frame, dt, eq_time, title, xlabel, ylabel, labels, linestyles, colors, plot_avg, last_frame, array_filename, gamma_A, gamma_B, phi_A, phi_B, N_A, N_B):

		#Analysis properties
		self.init_frame = init_frame
		self.final_frame = final_frame
		self.dt = dt
		self.eq_time = int(eq_time / dt)

		self.title = title
		self.xlabel = xlabel
		self.ylabel = ylabel
		self.labels = labels.split()

		self.linestyles = linestyles.split()
		self.colors = colors.split()

		self.plot_avg = plot_avg
		self.last_frame = last_frame

		self.array_filename = array_filename

		#System properties
		self.gamma_A = gamma_A
		self.gamma_B = gamma_B
		self.gamma_AB = np.sqrt(self.gamma_A * self.gamma_B)

		self.phi_A = phi_A
		self.phi_B = phi_B

		self.N_A = N_A
		self.N_B = N_B

		self.q_theo = 2 * np.sqrt(self.phi_A) * np.exp(self.gamma_AB - 0.5 * self.gamma_A - 0.5)

	#Analyse data
	def aggregate_analysis(self):

	    aggregates_T = np.load(self.array_filename)

	    N_aggregates_t = []
	    Average_lenght_aggregates_t = []

	    os.chdir(cluster_folder)
	    f = open('aggregate_analysis.txt', 'w') #Save data in file for more customizable plots
	    f.write('t\tN_aggregates\tavg_length\n')

	    if self.last_frame == True:
	    	self.final_frame = len(aggregates_T)

	    if self.final_frame > len(aggregates_T):
	    	self.final_frame = len(aggregates_T)

	    for i in range(self.init_frame, self.final_frame, self.dt):
	        N_aggregates_t.append(len(aggregates_T[i]))
	        Average_lenght_aggregates_t.append(np.mean(aggregates_T[i]))

	        f.write(str(i) + '\t' + str(len(aggregates_T[i])) + '\t' + str(round(np.mean(aggregates_T[i]), 2)) + '\n')

	    #Average number of aggregates in time (from equilibrium)
	    avg_N = np.mean(N_aggregates_t[self.eq_time :])

	    avg_N_aggregates = np.linspace(avg_N, avg_N, len(N_aggregates_t))

	    #Average average length of aggregates in time
	    avg_l = np.mean(Average_lenght_aggregates_t[self.eq_time :])

	    avg_avg_length = np.linspace(avg_l, avg_l, len(Average_lenght_aggregates_t))

	    f.write('avg_N=%.2f\tavg_l=%.2f' % (avg_N, avg_l))
	    f.close()

	    #Time linspace (x axis)
	    t = np.linspace(0, len(aggregates_T), len(N_aggregates_t))

	    return t, N_aggregates_t, avg_N_aggregates, avg_N, Average_lenght_aggregates_t, avg_avg_length, avg_l, self.final_frame

	def linear_chains_statistics(self):
	        
	    N_B_in_linear_chains = np.load(self.array_filename)

	    B_in_chains = []
	    qs = []

	    os.chdir(linear_chains_folder)
	    f = open('linear_chains_analysis.txt', 'w') #Save data in file for more customizable plots
	    f.write('t\tN_B_linear_chains\tqs\n')

	    if self.last_frame == True:
	    	self.final_frame = len(N_B_in_linear_chains)

	    if self.final_frame > len(N_B_in_linear_chains):
	    	self.final_frame = len(N_B_in_linear_chains)

	    for i in range(self.init_frame, self.final_frame, self.dt):
	    	B_in_chains.append(N_B_in_linear_chains[i])

	    	q = N_B_in_linear_chains[i] / (self.N_B - N_B_in_linear_chains[i])

	    	qs.append(q)

	    	f.write(str(i) + '\t' + str(N_B_in_linear_chains[i]) + '\t' + str(q) +  '\n')

	    t = np.linspace(0, len(N_B_in_linear_chains), len(B_in_chains))

	    avg_B_in_chains = np.mean(B_in_chains[self.eq_time :])
	    avg_q = np.mean(qs[self.eq_time :])

	    f.write('avg_B_in_chains=%.2f\tavg_q=%.2f' % (avg_B_in_chains, avg_q))
	    f.close()

	    return t, B_in_chains, avg_B_in_chains, qs, avg_q, self.final_frame

	def lateral_chains_statistics(self, filename1, filename2):
		N_B_in_linear_chains_1 = np.load(filename1)
		N_B_in_linear_chains_2 = np.load(filename2)

		considered_length = len(N_B_in_linear_chains_1)
		different_length = False

		if len(N_B_in_linear_chains_1) != len(N_B_in_linear_chains_2):
			considered_length = min(len(N_B_in_linear_chains_1), len(N_B_in_linear_chains_2))
			different_length = True

		B_in_chains_1 = []
		B_in_chains_2 = []
		N_lateral = []
		qs = []

		if self.last_frame == True:
			self.final_frame = considered_length

		if self.final_frame > considered_length:
			self.final_frame = considered_length

		f = open('lateral_analysis.txt', 'w')
		f.write('t\tN_B_1\tN_B_1\tN_lateral\tq_lat_lin\n')

		for i in range(self.init_frame, self.final_frame, self.dt):
			B_in_chains_1.append(N_B_in_linear_chains_1[i])
			B_in_chains_2.append(N_B_in_linear_chains_2[i])

			N_lat = N_B_in_linear_chains_2[i] - N_B_in_linear_chains_1[i]
			N_lateral.append(N_lat)

			q_lat_lin = N_B_in_linear_chains_2[i] / (self.N_B - N_B_in_linear_chains_2[i]) #q due lateral and linear chains
			qs.append(q_lat_lin)


			f.write(str(i) + '\t' + str(N_B_in_linear_chains_1[i]) + '\t' + str(N_B_in_linear_chains_2[i]) + '\t' + str(N_lat) + '\t' + str(q_lat_lin) +  '\n')

		t = np.linspace(0, len(N_B_in_linear_chains_1), len(B_in_chains_1))

		avg_B_in_chains_1 = np.mean(B_in_chains_1[self.eq_time :])
		avg_B_in_chains_2 = np.mean(B_in_chains_2[self.eq_time :])
		avg_N_lateral = np.mean(N_lateral[self.eq_time :])
		avg_q = np.mean(qs[self.eq_time :])

		f.write('avg_B_in_chains_1=%.2f\tavg_B_in_chains_2=%.2f\tavg_N_lateral=%.2f\tavg_q=%.2f' % (avg_B_in_chains_1, avg_B_in_chains_2, avg_N_lateral, avg_q))
		f.close()

		return t, N_lateral, avg_N_lateral, qs, avg_q, self.final_frame, different_length

	def general_aggregates_analysis(self):
		general_aggregates = np.load(self.array_filename)

		if self.last_frame == True:
			self.final_frame = len(general_aggregates)

		if self.final_frame > len(general_aggregates):
			self.final_frame = len(general_aggregates)

		general_aggregates_new = []

		os.chdir(general_aggregates_folder)
		f = open('general_aggregate_analysis.txt', 'w')
		f.write('t\tN_aggregates\n')

		for i in range(self.init_frame, self.final_frame, self.dt):
			general_aggregates_new.append(general_aggregates[i])	    
			
			f.write(str(i) + '\t' + str(general_aggregates[i]) +  '\n')

		t = np.linspace(self.init_frame, self.final_frame, len(general_aggregates_new))
		avg = np.mean(general_aggregates_new[self.eq_time:])

		f.write('avg=%.2f' % avg)
		f.close()

		return t, general_aggregates_new, avg, self.final_frame

	def isolated_particles_analysis(self):
		isolated_particles = np.load(self.array_filename)

		if self.last_frame == True:
			self.final_frame = len(isolated_particles)

		if self.final_frame > len(isolated_particles):
			self.final_frame = len(isolated_particles)

		isolated_particles_new = []

		os.chdir(isolated_particles_folder)
		f = open('general_aggregate_analysis.txt', 'w')
		f.write('t\tisolated_part\n')

		for i in range(self.init_frame, self.final_frame, self.dt):
			isolated_particles_new.append(isolated_particles[i])	    

			f.write(str(i) + '\t' + str(isolated_particles[i]) +  '\n')

		t = np.linspace(self.init_frame, self.final_frame, len(isolated_particles_new))
		avg = np.mean(isolated_particles_new[self.eq_time:])

		f.write('avg=%.2f' % avg)
		f.close()

		return t, isolated_particles_new, avg, self.final_frame

	def isolated_aggregates_analysis(self):

		isolated_aggregates = np.load(self.array_filename)

		if self.last_frame == True:
			self.final_frame = len(isolated_aggregates)

		if self.final_frame > len(isolated_aggregates):
			self.final_frame = len(isolated_aggregates)

		isolated_aggregates_new = []

		os.chdir(isolated_aggregates_folder)
		f = open('isolated_aggregates_analysis.txt', 'w')
		f.write('t\tB_in_isolated_aggregates\n')

		for i in range(self.init_frame, self.final_frame, self.dt):
			isolated_aggregates_new.append(isolated_aggregates[i])	    

			f.write(str(i) + '\t' + str(isolated_aggregates[i]) +  '\n')

		t = np.linspace(self.init_frame, self.final_frame, len(isolated_aggregates_new))
		avg = np.mean(isolated_aggregates_new[self.eq_time:])

		f.write('avg=%.2f' % avg)
		f.close()

		return t, isolated_aggregates_new, avg, self.final_frame


	#Plot analysed data
	def N_clusters(self, t, N_aggregates_t, avg_N_aggregates, avg_N, filename, default):

		if default == True: #Set default settings
			plt.plot(t, N_aggregates_t)

			if self.plot_avg == True: 
				plt.plot(t, avg_N_aggregates)

		else:
			plt.plot(t, N_aggregates_t, ls=self.linestyles[0], color=self.colors [0], label=self.labels[0])

			if self.plot_avg == True: 
				plt.plot(t, avg_N_aggregates, ls=self.linestyles[1], color=self.colors[1], label=self.labels[1] + ' ' + str(int(round(avg_N, 0))))

		plt.title(self.title)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)

		plt.legend()
		plt.show()
		plt.gcf().clear()
		plt.close('all')

	def avg_cluster_length(self, t, Average_lenght_aggregates_t, avg_avg_length, avg_l, filename, default):

		if default == True: #Set default settings
			plt.plot(t, Average_lenght_aggregates_t)

			if self.plot_avg == True:
				plt.plot(t, avg_avg_length)

		else:

			plt.plot(t, Average_lenght_aggregates_t, ls=self.linestyles[0] ,color=self.colors [0], label=self.labels[0])

			if self.plot_avg == True:
				plt.plot(t, avg_avg_length, ls=self.linestyles[1], color=self.colors[1], label=self.labels[1] + ' ' + str(round(avg_l, 2)))

		plt.title(self.title)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)

		plt.legend()
		plt.show()
		plt.gcf().clear()
		plt.close('all')

	def N_B_in_linear_chains(self, t, N_B_in_chains, avg_B_in_chains, filename, default):
		
		if default == True:
			plt.plot(t, N_B_in_chains)

			if self.plot_avg == True:
				plt.plot(t, np.linspace(avg_B_in_chains, avg_B_in_chains, len(N_B_in_chains)))

		else:
			plt.plot(t, N_B_in_chains, ls=self.linestyles[0], color=self.colors[0], label=self.labels[0])

			if self.plot_avg == True:
				plt.plot(t, np.linspace(avg_B_in_chains, avg_B_in_chains, len(N_B_in_chains)), ls=self.linestyles[1], color=self.colors[1], label=self.labels[1] + ' ' + str(int(round(avg_B_in_chains))))

		plt.title(self.title)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)

		plt.legend()
		plt.show()
		plt.gcf().clear()
		plt.close('all')

	def q_linear_chains(self, t, qs, avg_q, filename, default):

		if default == True:

			plt.plot(t, qs)
			plt.plot(t, np.linspace(self.q_theo, self.q_theo, len(qs)), ls='-', color='g', label='Theoretical q value')

			if self.plot_avg == True:
				plt.plot(t, np.linspace(avg_q, avg_q, len(qs)))

		else:
			plt.plot(t, qs, ls=self.linestyles[0], color=self.colors[0], label=self.labels[0])
			plt.plot(t, np.linspace(self.q_theo, self.q_theo, len(qs)), ls='-', color='g', label='Theoretical q value %.3f' % self.q_theo)

			if self.plot_avg == True:
				plt.plot(t, np.linspace(avg_q, avg_q, len(qs)), ls=self.linestyles[1], color=self.colors[1], label=self.labels[1] + ' ' + str(round(avg_q, 3)))

		plt.title(self.title)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)

		plt.legend()
		plt.show()
		plt.gcf().clear()
		plt.close('all')

	def N_general_aggregates(self, t, N_aggregates, avg, filename, default):

		if default == True:
			plt.plot(t, N_aggregates)

			if self.plot_avg == True:
				plt.plot(t, np.linspace(avg, avg, len(N_aggregates)))

		else:
			plt.plot(t, N_aggregates, ls=self.linestyles[0], color=self.colors[0], label=self.labels[0])

			if self.plot_avg == True:
				plt.plot(t, np.linspace(avg, avg, len(N_aggregates)), ls=self.linestyles[1], color=self.colors[1], label=self.labels[1] + ' ' + str(int(round(avg))))

		plt.title(self.title)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)

		plt.legend()
		plt.show()
		plt.gcf().clear()
		plt.close('all')

	def N_isolated_particles(self, t, isolated_particles, avg, filename, default):

		if default == True:
			plt.plot(t, isolated_particles)

			if self.plot_avg == True:
				plt.plot(t, np.linspace(avg, avg, len(isolated_particles)))

		else:
			plt.plot(t, isolated_particles, ls=self.linestyles[0], color=self.colors[0], label=self.labels[0])

			if self.plot_avg == True:
				plt.plot(t, np.linspace(avg, avg, len(isolated_particles)), ls=self.linestyles[1], color=self.colors[1], label=self.labels[1] + ' ' + str(int(round(avg))))

		plt.title(self.title)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)

		plt.legend()
		plt.show()
		plt.gcf().clear()
		plt.close('all')

	def plot_rdf(self, filename_load, filename_save, default):

		rdf_data = np.load(filename_load)

		x = rdf_data[0]
		y = rdf_data[1]

		if default == True:
			plt.plot(x, y)

			if self.plot_avg == True:
				plt.plot(x, np.linspace(1, 1, len(x)))

		else:
			plt.plot(x, y, ls=self.linestyles[0], color=self.colors[0], label=self.labels[0])

			if self.plot_avg == True:
				plt.plot(x, np.linspace(1, 1, len(x)), ls=self.linestyles[1], color=self.colors[1], label=self.labels[1])

		plt.title(self.title)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)

		plt.legend()
		plt.show()
		plt.gcf().clear()
		plt.close('all')


app = QApplication(sys.argv)
_window=Window()
_window.show()
app.exec_()
