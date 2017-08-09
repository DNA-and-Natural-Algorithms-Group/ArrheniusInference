from __future__ import division, absolute_import
import ConfigParser
import numpy as np
import statistics
from scipy import stats
from scipy.optimize import minimize
import math
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import cPickle as pickle
import random
import multiprocessing
from multiprocessing import Lock, Process, Queue, current_process, Manager
import timeit
import myenums
__all__ = ['minimize', 'minimize_scalar']
import learndnakinetics


from matplotlib.ticker import ScalarFormatter
from matplotlib.font_manager import FontProperties
import os
from collections import OrderedDict

MEGA = 1000000
ORDER_OF_MAGNITUDE = 0.4771
USE_ONLY_FINALSTEP = True #When use_only_finalstep  = True , only  the last step of each walker is used!


class ForMultiProcess(object):
	def __init__(self, functionName, arguments):
		self.functionName = functionName
		self.arguments = arguments



class PlotOptions(object):
	"""Setting options for different plots based on the plots in the literature. """
	def __init__(self, plot_errorbars, use_ylim, use_yticks, logscale, fig_name, document, temps_name_type, temps_name,
				 experimental_property_type, experimental_property_cell, type, colors, style, _zip, xticks, yticks, ylim,
				 xlabel, ylabel, title, loc, markersize, linewidth, title_fontsize, label_fontsize, propsize, capsize,
				 show_legend, labelsize):
		self.use_yticks = use_yticks
		self.plot_errorbars = plot_errorbars
		self.use_ylim = use_ylim
		self.logscale = logscale
		self.document = document
		self.temps_name_type = temps_name_type
		self.temps_name = temps_name
		self.experimental_property_cell = experimental_property_cell
		self.experimental_property_type = experimental_property_type
		self.type = type
		self.style = style
		self.colors = colors
		self.fig_name = fig_name
		self._zip = _zip
		self.xticks = xticks
		self.yticks = yticks
		self.ylim = ylim
		self.xlabel = xlabel
		self.ylabel = ylabel
		self.loc = loc
		self.title = title
		self.label_fontsize = label_fontsize
		self.markersize = markersize
		self.linewidth = linewidth
		self.title_fontsize = title_fontsize
		self.linewidth = linewidth
		self.propsize = propsize
		self.capsize = capsize
		self.show_legend = show_legend
		self.labelsize = labelsize

def find_mean_errorbars(predictions_ensemble):
	"""returns the mean of the predictions  in predictions_ensemble. errorbar_min and errorbar_up  indicate the errorbar """
	predicted_property = average = np.average(predictions_ensemble)
	errorbar_min = average - min(predictions_ensemble)
	errorbar_up = max(predictions_ensemble) - average
	return predicted_property, (errorbar_min, errorbar_up )


def find(theta, counter_cell, document, row, experimental_property_type, experimental_property_cell, type):

	#  Plots of different literature represent different properties. Therefore there are different options in this function.
	predictions_ensemble = []
	predicted_logreactionrateconstant = []
	
	for th in theta:
		th = str(th)
		plrc = predicted_logreactionrateconstants[th, counter_cell, document]
		experimental_logreactionrateconstant= experimental_logreactionrateconstants[th, counter_cell, document]
		keffprob = plot_probs[th]
		if keffprob != -np.inf:
			if type == 1:
				predictions_ensemble.append(plrc)
			elif type == 2:

				predictions_ensemble.append(10 ** plrc)
			elif type == 3:
				predictions_ensemble.append(math.pow(10.0, plrc) / MEGA)
			elif type == 4:
				predictions_ensemble.append((1 / math.pow(10, plrc)) * MEGA)

			predicted_logreactionrateconstant.append(plrc)
	predicted_property, variance = find_mean_errorbars(predictions_ensemble)
	
	
	if experimental_property_type == 1:
		experimental_property = float(row[counter_cell][experimental_property_cell])
	elif experimental_property_type == 2:
		experimental_property = float(row[counter_cell][experimental_property_cell]) / MEGA
	elif experimental_property_type == 3:
		experimental_property = float(row[counter_cell][experimental_property_cell]) * MEGA
	elif experimental_property_type == 4:
		experimental_property = float(row[counter_cell][experimental_property_cell])
		experimental_property = np.log(experimental_property) / np.log(10)
	return experimental_property, predicted_property, variance, experimental_logreactionrateconstant, predicted_logreactionrateconstant


def obtain_points(theta, _zip, document, temps_name_type, temps_name, experimental_property_type, experimental_property_cell, type):
	
	success_train = 0
	success_test = 0
	n_tests = 0
	n_trains = 0
	infis = 0
	names = []
	temps = dict()
	predicted_properties = dict()
	experimental_properties = dict()
	errorbars = dict()
	errorbars1 = dict()
	errorbars2 = dict()

	row = learndnakinetics.open_document(document)
	test_error = 0
	train_error = 0
	try:
		
		for counter_cell in range(1, learndnakinetics.n_csv_rows(row)):
			name = str(row[counter_cell][0].rstrip())
			if name not in names:
				temps[name] = []
				predicted_properties[name] = []
				experimental_properties[name] = []
				errorbars1[name] = []
				errorbars2[name] = []
				errorbars[name] = []
				names.append(name)

			if temps_name_type == 0:
				mo = int(row[counter_cell][temps_name])
				if mo == -1:
					mo = 1
				temps[name].append(mo)

			if temps_name_type == 1:
				temps[name].append(float(row[counter_cell][temps_name]))

			elif temps_name_type == 2:
				temps[name].append(len(row[counter_cell][temps_name].rstrip()))

			elif temps_name_type == 3:
				temps[name].append(counter_cell)

			elif temps_name_type == 4:
				temps[name].append((1000 / (float(row[counter_cell][temps_name]) + 273.15)))

			experimental_property, predicted_property, variance, experimental_logreactionrateconstant, predicted_logreactionrateconstant = find(
				theta, counter_cell, document, row, experimental_property_type, experimental_property_cell, type)

			predicted_properties[name].append(predicted_property)
			experimental_properties[name].append(experimental_property)
			errorbars1[name].append(variance[0])
			errorbars2[name].append(variance[1])

			error = abs(np.average(predicted_logreactionrateconstant) - experimental_logreactionrateconstant) ** 2
			if error <= ORDER_OF_MAGNITUDE and counter_cell in traintestset[document, myenums.SetName.TRAIN.value]:
				success_train += 1
			if error <= ORDER_OF_MAGNITUDE and counter_cell in traintestset[document, myenums.SetName.TEST.value]:
				if '/hairpin4/Fig5' in document and counter_cell == 5:
					print "Not including temperature 20 degree celcius because it already  exists in the training set"
				else:
					success_test += 1

			if counter_cell in traintestset[document, myenums.SetName.TEST.value]:
				if error != np.inf:
					if '/hairpin4/Fig5' in document and counter_cell == 5:
						print " not including temperature 20 degree celcius because already in training "
					else:
						test_error += error
				else:
					infis += 1
				if '/hairpin4/Fig5' in document and counter_cell == 5:
					print " not including temperature 20 degree for /hairpin4/fig5  because  it already  exists in the training set"
				else:
					n_tests += 1
			else:
				if error != np.inf:
					train_error += error
				else:
					infis += 1
				n_trains += 1

	except ValueError:
		raise ValueError('Exception:(!')
	for name in names:
		errorbars[name] = (errorbars1[name], errorbars2[name])

	aa = np.zeros((7))
	aa[0] = train_error
	aa[1] = success_train
	aa[2] = n_trains
	aa[3] = test_error
	aa[4] = success_test
	aa[5] = n_tests
	aa[6] = infis
	done_queue = (aa)
	if n_tests != 0:
		test_error_mean = test_error / n_tests
	else:
		test_error_mean = "No Test"
	if n_trains != 0:
		train_error_mean = test_error / n_trains
	else:
		train_error_mean = "No Train"
	return [names, temps, predicted_properties, test_error_mean, train_error_mean, experimental_properties, errorbars,
			done_queue]


def draw_plot(done_queue_results, theta, po):
	#draws plots
	done = 0
	fig = plt.figure()
	errors_closeopen_test = dict()
	errors_closeopen_train = dict()

	for _zip in po._zip:

		if '/hairpin4/Fig' in po.document:

			if _zip == False:
				po.style = ['s', 's', 's']
			else:
				po.style = ['o', 'o', 'o']

		if '.csv' in po.document:
			document = po.document
		else:
			document = po.document + str(int(_zip)) + '.csv'

		[names, temps, predicted_properties, test_error, train_error, experimental_properties, errorbars,
		 done_queue] = obtain_points(theta, _zip, document, po.temps_name_type, po.temps_name, po.experimental_property_type,
									po.experimental_property_cell, po.type)
		i = 0
		done += done_queue
		errors_closeopen_test[_zip] = test_error
		errors_closeopen_train[_zip] = train_error

		for obj in names:
			color = po.colors[i]
			prediction_linestyle = '--'
			real_linestyle = '-'
			if _zip == True:
				if po.plot_errorbars == True:
					plt.errorbar(temps[obj], predicted_properties[obj], yerr=[errorbars[obj][0], errorbars[obj][1]],
								 capsize=po.capsize, linestyle=prediction_linestyle, marker=po.style[i],
								 color=color, linewidth=po.linewidth, markersize=po.markersize)
				else:
					plt.plot(temps[obj], predicted_properties[obj], linestyle=prediction_linestyle, marker=po.style[i],
							 color=color, linewidth=po.linewidth, markersize=po.markersize)

				plt.plot(temps[obj], experimental_properties[obj], linestyle=real_linestyle, marker=po.style[i],
						 color=po.colors[i], linewidth=po.linewidth, label=obj, markersize=po.markersize)
			else:

				if po.plot_errorbars == True:
					plt.errorbar(temps[obj], predicted_properties[obj], capsize=po.capsize, linestyle=':',
								 yerr=[errorbars[obj][0], errorbars[obj][1]], marker=po.style[i], color=po.colors[i],
								 markerfacecolor='none', markersize=po.markersize)
				else:
					plt.plot(temps[obj], predicted_properties[obj], linestyle=':', marker=po.style[i],
							 color=po.colors[i], markerfacecolor='none', markersize=po.markersize)

				if len(po._zip) == 1:
					plt.plot(temps[obj], experimental_properties[obj], linestyle='-', marker=po.style[i],
							 color=po.colors[i], label=obj, markerfacecolor='none', markersize=po.markersize)
				else:
					plt.plot(temps[obj], experimental_properties[obj], linestyle='-', marker=po.style[i],
							 color=po.colors[i], markerfacecolor='none', markersize=po.markersize)

			i += 1

	if po.logscale == True:
		plt.yscale('log')
	plt.xticks(po.xticks)
	plt.tick_params(axis='both', which='major', labelsize=po.labelsize)
	if learndnakinetics.DATASET_PATH + '/helix3/zhang_1.csv' in document:
		fig = matplotlib.pyplot.gcf()
		fig.set_size_inches(30, 10.5)

	if po.use_ylim == True:
		axes = plt.gca()
		axes.set_ylim(po.ylim)
	if po.use_yticks == True:
		plt.yticks(po.yticks)
	# ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	plt.title(po.title, fontsize=po.title_fontsize)
	if po.show_legend == True:
		plt.legend(loc=po.loc, borderaxespad=0., prop={'size': po.propsize})
	plt.xlabel(po.xlabel, fontsize=po.label_fontsize)
	plt.ylabel(po.ylabel, fontsize=po.label_fontsize)
	plt.savefig(po.fig_name, bbox_inches='tight')
	plt.show()
	plt.close(fig)
	done_queue_results.put(done)


def draw_plot_Machinek(done_queue_results, theta, po):
	#draws Fig. 2d from paper: "Programmable energy landscapes for kinetic control of DNA strand displacement".
	done = 0

	fig, ax = plt.subplots()
	errors_closeopen_test = dict()
	errors_closeopen_train = dict()

	[names, temps, predicted_properties, test_error, train_error, experimental_properties, errorbars,
	 done_queue] = obtain_points(theta, po._zip, po.document, po.temps_name_type, po.temps_name,
								po.experimental_property_type, po.experimental_property_cell, po.type)
	done += done_queue
	i = 0
	errors_closeopen_test[po._zip] = test_error
	errors_closeopen_train[po._zip] = train_error
	for obj in names:

		if po._zip == True:

			ak = predicted_properties[obj][1:]
			at = temps[obj][1:]
			if po.plot_errorbars == True:

				plt.errorbar(at, ak, yerr=[errorbars[obj][0][1:], errorbars[obj][1][1:]], linestyle='--',
							 marker=po.style[i], color=po.colors[i], linewidth=po.linewidth, markersize=po.markersize,
							 capsize=po.capsize)
			else:
				plt.plot(at, ak, linestyle='--', marker=po.style[i], color=po.colors[i], linewidth=po.linewidth,
						 markersize=po.markersize)
			ax.set_yscale("log", nonposy='clip')
			ax.annotate('', xy=(1, predicted_properties[obj][0]), xytext=(2, predicted_properties[obj][0]),
						arrowprops=dict(arrowstyle="-|>", ls='dashed', color=po.colors[i]), )
			plt.plot(temps[obj][1:], experimental_properties[obj][1:], linestyle='-', marker=po.style[i],
					 color=po.colors[i], linewidth=po.linewidth, label=obj, markersize=po.markersize)
			ax.set_yscale("log", nonposy='clip')
			ax.annotate('', xy=(1, experimental_properties[obj][0]), xytext=(2, experimental_properties[obj][0]),
						arrowprops=dict(arrowstyle="-|>", color=po.colors[i]), )

		i += 1

	plt.title(po.title, fontsize=po.title_fontsize)

	ax.xaxis.set_major_formatter(ScalarFormatter())
	plt.xticks(po.xticks)
	xticks = ax.xaxis.get_major_ticks()
	xticks[0].label1.set_visible(False)
	plt.tick_params(axis='both', which='major', labelsize=po.labelsize)

	axes = plt.gca()
	axes.set_ylim(po.yticks)

	if po.show_legend == True:
		plt.legend(loc=po.loc, borderaxespad=0., prop={'size': po.propsize})
	plt.xlabel(po.xlabel, fontsize=po.label_fontsize)
	plt.ylabel(po.ylabel, fontsize=po.label_fontsize)
	plt.savefig(po.fig_name, bbox_inches='tight')
	plt.show()
	plt.close(fig)
	
	done_queue_results.put(done)


def main():
	global  traintestset, predicted_logreactionrateconstants, experimental_logreactionrateconstants, plot_probs
	learndnakinetics.set_configuration() #setting some configuations
	learndnakinetics.use_all_data = True
	configParser = ConfigParser.ConfigParser()
	configParser.readfp(open(r'config_file.txt'))
	CONFIG_NAME = 'plot'
	learned_parameterset = configParser.get(CONFIG_NAME, 'learned_parameterset')
	plot_errorbars = False
	if  learned_parameterset==   myenums.LearnedParameters.ARRHENIUSMCMC.value[0]  or learned_parameterset ==  myenums.LearnedParameters.METROPOLISMCMC.value[0]   :
		plot_errorbars =True #For the ensemble MCMC draw errorbars!
	n_processors = configParser.getint(CONFIG_NAME, 'n_processors_plot')
	if not os.path.exists("plot/" + str(learned_parameterset)):
		os.makedirs("plot/" + str(learned_parameterset))

	traintestset, predicted_logreactionrateconstants, experimental_logreactionrateconstants, plot_probs, title_name = learndnakinetics.plot_rates(configParser.get(CONFIG_NAME, 'pkl_file_plot'),  learned_parameterset , USE_ONLY_FINALSTEP)
	theta = []
 
	for t in plot_probs:
		theta.append(t)
	if learned_parameterset == myenums.LearnedParameters.METROPOLISMAP.value[0]:
		# According to our DNA23 paper, only show legends for the top left figure.
		show_legend = True
	else:
		show_legend = False

	markersize = 8
	linewidth = 2
	title_fontsize = 23
	label_fontsize = 19
	propsize = 17
	capsize = 4
	labelsize = 17

	done_queue_results = multiprocessing.Manager().Queue()
	dataset_list = []

	filename = "plot/" + str(learned_parameterset)
	if not os.path.exists(filename):
		os.makedirs(filename)


	po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=True, logscale=False,
					 fig_name=filename + '/machinek2014.pdf',
					 document=learndnakinetics.DATASET_PATH + '/three_waystranddisplacement2/Fig2.csv', temps_name_type=0, temps_name=2,
					 experimental_property_type=1, experimental_property_cell=9, type=2, colors=['black', 'r', 'g'],
					 style=['s', 'o', 'v'], _zip=True, xticks=[1, 2, 4, 6, 8, 10, 12, 14, 16],
					 yticks=[10 ** 2, 10 ** 8], ylim=[10 ** 2, 10 ** 8], xlabel='Mismatch position (nt)',
					 ylabel=r"$\mathregular{K(M^{-1}s^{-1})}$", loc=4, markersize=markersize, linewidth=linewidth,
					 title_fontsize=title_fontsize, label_fontsize=label_fontsize, propsize=propsize, capsize=capsize,
					 show_legend=show_legend, labelsize=labelsize)
	dataset_list.append(ForMultiProcess(draw_plot_Machinek , (done_queue_results, theta, po)))

	
	for j in ['a', 'b']:
		po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=True, logscale=True,
						 _zip=[True, False], fig_name=filename + '/kimfig5' + str(j) + '.pdf',
						 document=learndnakinetics.DATASET_PATH + '/hairpin4/Fig5' + str(j) + '_', temps_name_type=1, temps_name=3,
						 experimental_property_type=2, experimental_property_cell=5, type=3, colors=['k', 'k'], style="",
						 xticks=np.arange(3.1, 3.63, 0.1), yticks=[0.1, 1, 10], ylim=[0, 10],
						 xlabel=r"$\mathregular{1000/T (K^{-1})}$", ylabel=r"$\mathregular{k_{op/cl}(10^6s^{-1})}$",
						 loc=0, markersize=markersize, linewidth=linewidth, title_fontsize=title_fontsize,
						 label_fontsize=label_fontsize, propsize=propsize, capsize=capsize, show_legend=show_legend,
						 labelsize=labelsize)
		dataset_list.append(ForMultiProcess(draw_plot, (done_queue_results, theta, po)))
	
	for j in [4, 6]:
		po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=True, logscale=True,
						 _zip=[True, False], fig_name=filename + '/bonnet1998' + str(j) + '.pdf',
						 document=learndnakinetics.DATASET_PATH + '/hairpin/Fig' + str(j) + '_', temps_name_type=1, temps_name=3,
						 experimental_property_type=1, experimental_property_cell=5, type=2,
						 colors=['r', 'g', 'b', 'y', 'k'], style=['o', 's', 'D', '^', 'p', '1'],
						 xticks=np.arange(3.1, 3.63, 0.1), yticks=[10, 100, 1000, 10000, 100000],
						 ylim=[10, 2 * 10 ** 5], xlabel=r"$\mathregular{1000/T (K^{-1})}$",
						 ylabel=r"$\mathregular{k_{-},k_{+}(s^{-1})}$", loc=0, markersize=markersize,
						 linewidth=linewidth, title_fontsize=title_fontsize, label_fontsize=label_fontsize,
						 propsize=propsize, capsize=capsize, show_legend=show_legend, labelsize=labelsize)
		dataset_list.append(ForMultiProcess(draw_plot, (done_queue_results, theta, po)))

	po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=True, logscale=False,
					 _zip=[True, False], fig_name=filename + '/morrison.pdf', document=learndnakinetics.DATASET_PATH + '/helix/Fig6_',
					 temps_name_type=1, temps_name=3, experimental_property_type=1, experimental_property_cell=5, type=1,
					 colors=['r', 'g', 'b', 'm', 'y', 'k'], style=['*', 'o', 'D', 's', 'x', '^', '*'],
					 xticks=np.arange(2.9, 3.7, 0.2), yticks=np.arange(-3, 8, 2), ylim=[-4, 8],
					 xlabel=r"$\mathregular{1000/T (K^{-1})}$", ylabel=r"$\mathregular{log(K)}$", loc=0,
					 markersize=markersize, linewidth=linewidth, title_fontsize=title_fontsize,
					 label_fontsize=label_fontsize, propsize=propsize, capsize=capsize, show_legend=show_legend,
					 labelsize=labelsize)
	dataset_list.append(ForMultiProcess(draw_plot, (done_queue_results, theta, po)))

	po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=True, logscale=False,
					 _zip=[False], fig_name=filename + '/ReynaldoDisassocation.pdf',
					 document=learndnakinetics.DATASET_PATH + '/helix1/Fig6a.csv', temps_name_type=4, temps_name=3,
					 experimental_property_type=4, experimental_property_cell=5, type=1,
					 colors=['r', 'g', 'b', 'm', 'y', 'k'], style=['*', 'o', 'D', 's', 'x', '^', '*'],
					 xticks=np.arange(3, 3.35, 0.05), yticks=np.arange(-5.5, -2, 0.5), ylim=[-6.5, -2],
					 xlabel=r"$\mathregular{1000/T (K^{-1})}$", ylabel=r"$\mathregular{log(K)}$", loc=0,
					 markersize=markersize, linewidth=linewidth, title_fontsize=title_fontsize,
					 label_fontsize=label_fontsize, propsize=propsize, capsize=capsize, show_legend=show_legend,
					 labelsize=labelsize)
	dataset_list.append(ForMultiProcess(draw_plot, (done_queue_results, theta, po)))

	po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=True, logscale=False,
					 _zip=[True], fig_name=filename + '/ReynaldoSequential.pdf',
					 document=learndnakinetics.DATASET_PATH + '/three_waystranddisplacement1/Fig6b.csv', temps_name_type=4, temps_name=3,
					 experimental_property_type=4, experimental_property_cell=5, type=1,
					 colors=['r', 'g', 'b', 'm', 'y', 'k'], style=['*', 'o', 'D', 's', 'x', '^', '*'],
					 xticks=np.arange(3, 3.35, 0.05), ylim=[0, 3], yticks=np.arange(0, 3, 0.5),
					 xlabel=r"$\mathregular{1000/T (K^{-1})}$", ylabel=r"$\mathregular{log(K)}$", loc=0,
					 markersize=markersize, linewidth=linewidth, title_fontsize=title_fontsize,
					 label_fontsize=label_fontsize, propsize=propsize, capsize=capsize, show_legend=show_legend,
					 labelsize=labelsize)
	dataset_list.append(ForMultiProcess(draw_plot, (done_queue_results, theta, po)))

	po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=False, logscale=False,
					 _zip=[True], fig_name=filename + '/4wayNadine.pdf',
					 document=learndnakinetics.DATASET_PATH + '/four_waystrandexchange/Table5.2.csv', temps_name_type=1, temps_name=2,
					 experimental_property_type=4, experimental_property_cell=8, type=1,
					 colors=['r', 'g', 'b', 'm', 'y', 'k'], style=['*', 'o', 'D', 's', 'x', '^', '*'],
					 xticks=np.arange(0, 9, 2), yticks="", ylim=[-3, 8], xlabel="n (nt)", ylabel="log(K)", loc=0,
					 markersize=markersize, linewidth=linewidth, title_fontsize=title_fontsize,
					 label_fontsize=label_fontsize, propsize=propsize, capsize=capsize, show_legend=show_legend,
					 labelsize=labelsize)
	dataset_list.append(ForMultiProcess(draw_plot, (done_queue_results, theta, po)))

	po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=True, logscale=False,
					 fig_name=filename + '/zhang2009.pdf', document=learndnakinetics.DATASET_PATH + '/three_waystranddisplacement/Fig3b.csv',
					 temps_name_type=1, temps_name=1, experimental_property_type=1, experimental_property_cell=7, type=1,
					 colors=['black', 'r', 'g', 'y', 'b', 'p'], style=['*', 'o', 'D'], _zip=[True],
					 xticks=np.arange(0, 16, 5), yticks=np.arange(0, 10, 2), ylim=[0, 8], xlabel='Toehold length (nt)',
					 ylabel='log(K)', loc=0, markersize=markersize, linewidth=linewidth, title_fontsize=title_fontsize,
					 label_fontsize=label_fontsize, propsize=propsize, capsize=capsize, show_legend=show_legend,
					 labelsize=labelsize)
	dataset_list.append(ForMultiProcess(draw_plot, (done_queue_results, theta, po)))

	po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=True, logscale=True,
					 _zip=[True], fig_name=filename + '/kimTrue.pdf', document=learndnakinetics.DATASET_PATH + '/hairpin4/Table1_1.csv',
					 temps_name_type=2, temps_name=1, experimental_property_type=2, experimental_property_cell=5, type=3,
					 colors=['g', 'b', 'y', 'm', 'k'], style=['s', 'o', '^', 'v', 'D'], xticks=[0, 1, 2, 3],
					 yticks=[0.1, 1, 10], ylim=[0, 10], xlabel="# of dC-dG pairs", ylabel=r"$k_{cl}10^6s^{-1}$", loc=0,
					 markersize=markersize, linewidth=linewidth, title_fontsize=title_fontsize,
					 label_fontsize=label_fontsize, propsize=propsize, capsize=capsize, show_legend=show_legend,
					 labelsize=labelsize)
	dataset_list.append(ForMultiProcess(draw_plot, (done_queue_results, theta, po)))

	po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=True, logscale=True,
					 _zip=[False], fig_name=filename + '/kimFalse.pdf', document=learndnakinetics.DATASET_PATH + '/hairpin4/Table1_0.csv',
					 temps_name_type=2, temps_name=1, experimental_property_type=2, experimental_property_cell=5, type=3,
					 colors=['g', 'b', 'y', 'm', 'k'], style=['s', 'o', '^', 'v', 'D'], xticks=[0, 1, 2, 3],
					 yticks=[0.1, 1, 10], ylim=[0, 10], xlabel="# of dC-dG pairs", ylabel=r"$k_{op}10^6s^{-1}$", loc=0,
					 markersize=markersize, linewidth=linewidth, title_fontsize=title_fontsize,
					 label_fontsize=label_fontsize, propsize=propsize, capsize=capsize, show_legend=show_legend,
					 labelsize=labelsize)
	dataset_list.append(ForMultiProcess(draw_plot, (done_queue_results, theta, po)))

	po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=True, logscale=True,
					 _zip=[True], fig_name=filename + '/goddard2000T1.pdf',
					 document=learndnakinetics.DATASET_PATH + '/hairpin1/Fig3_T_1.csv', temps_name_type=1, temps_name=3,
					 experimental_property_type=3, experimental_property_cell=4, type=4, colors=['b', 'y', 'm', 'k'],
					 style=['x', '^', '*', 'D'], xticks=np.arange(3.1, 3.6, 0.1), yticks=[10, 100], ylim=[8, 10 ** 3],
					 xlabel='$1000/T (K^{-1})$', ylabel=r"$\mathregular{\tau_+  (\mu s)}$", loc=0,
					 markersize=markersize, linewidth=linewidth, title_fontsize=title_fontsize,
					 label_fontsize=label_fontsize, propsize=propsize, capsize=capsize, show_legend=show_legend,
					 labelsize=labelsize)
	dataset_list.append(ForMultiProcess(draw_plot, (done_queue_results, theta, po)))

	po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=True, logscale=True,
					 _zip=[False], fig_name=filename + '/goddard2000T0.pdf',
					 document=learndnakinetics.DATASET_PATH + '/hairpin1/Fig3_T_0.csv', temps_name_type=1, temps_name=3,
					 experimental_property_type=3, experimental_property_cell=4, type=4, colors=['b', 'y', 'm', 'k'],
					 style=['x', '^', '*', 'D'], xticks=np.arange(3.1, 3.50, 0.05), yticks=[10, 100, 1000, 10000],
					 ylim=[10, 10000], xlabel='$1000/T (K^{-1})$', ylabel=r"$\mathregular{\tau_-  (\mu s)}$", loc=0,
					 markersize=markersize, linewidth=linewidth, title_fontsize=title_fontsize,
					 label_fontsize=label_fontsize, propsize=propsize, capsize=capsize, show_legend=show_legend,
					 labelsize=labelsize)
	dataset_list.append(ForMultiProcess(draw_plot, (done_queue_results, theta, po)))

	po = PlotOptions(title=title_name, plot_errorbars=plot_errorbars, use_ylim=True, use_yticks=True, logscale=True,
					 _zip=[True], fig_name=filename + '/altanbonnet20034.pdf', document=learndnakinetics.DATASET_PATH + '/bubble/Fig4.csv',
					 temps_name_type=1, temps_name=4, experimental_property_type=3, experimental_property_cell=5, type=4,
					 colors=['b', 'r', 'g'], style=['^', 'o', 'v'], xticks=np.arange(3.1, 3.6, 0.1), yticks=[10, 100],
					 ylim=[0, 200], xlabel='$1000/T (K^{-1})$', ylabel=r"Timescale ($\mathregular{\mu s}$ )", loc=0,
					 markersize=markersize, linewidth=linewidth, title_fontsize=title_fontsize,
					 label_fontsize=label_fontsize, propsize=propsize, capsize=capsize, show_legend=show_legend,
					 labelsize=labelsize)
	dataset_list.append(ForMultiProcess(draw_plot, (done_queue_results, theta, po)))

	pool = multiprocessing.Pool(processes=n_processors)

	for ds in dataset_list:
		compile_error = pool.apply_async(ds.functionName, ds.arguments)
	#print compile_error.get()
	pool.close()
	pool.join()
	results = 0
	while not done_queue_results.empty():
		s = done_queue_results.get()
		print s
		results += s
	print results
	f = open('results' + str(learned_parameterset), 'w')
	f.write(str(results) + "\n")

if __name__ == "__main__":
	main()
