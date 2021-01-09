__author__ = 'lpeng'
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.basemap import Basemap, cm
import pandas as pd
import string

def CoherenceSinglePlot(coh, fs, freq, title):

	nvar = len(coh)
	col = ['Crimson', 'Orange', 'DeepSkyblue', 'SpringGreen', 'Purple']
	names = ['Tair','Rn','Wind','Humidity','VPD']


	[plt.semilogx(fs/freq, coh[i], color=col[i], linewidth=3, marker='o', markeredgecolor=col[i], alpha=0.3, label=names[i]) for i in range(nvar-1, -1, -1)]
	# plt.semilogx(fs/freq, coh, color=col[0], linewidth=2.5)
	plt.xlim([1, 2000])
	plt.text(1.5, 0.8, "%s" % title, fontsize=40)
	# plt.plot([365, 365], [0, 1], 'r--', linewidth=3.5)
	# plt.plot([120, 120], [0, 1], 'k--', linewidth=3.5)
	# plt.plot([30, 30], [0, 1], 'b--', linewidth=3.5)
	plt.grid(True)
	plt.tick_params(labelsize=26)
	plt.show()

	return

def CoherenceVarPlot(ax, coh, fs, freq):

	col = ['Crimson', 'Orange', 'DeepSkyblue', 'SpringGreen', 'Purple']
	# drivername = ['Tair','Rn','Wind','Humidity','VPD']
	drivername = [r'$T_a$', r'$R_n$', r'$u_2$',r'$e_a$',r'VPD']

	[ax.semilogx(fs/freq, coh[i, :], color=col[i], linewidth=2.5, alpha=0.6, label=drivername[i]) for i in (4, 1, 0, 3, 2)]  # xrange(nvar-1, -1, -1)

	## xlim, ylim, annotation
	ax.set_xlim([1, 2000]) # so that no space before and after the time series

	## add vertical lines
	ax.plot([365, 365], [0, 1], 'r--', linewidth=3)
	ax.plot([180, 180], [0, 1], 'k--', linewidth=3)
	ax.plot([90, 90], [0, 1], 'g--', linewidth=3)
	ax.plot([30, 30], [0, 1], 'b--', linewidth=3)
	ax.plot([7, 7], [0, 1], 'c--', linewidth=3)
	ax.plot([1, 2000], [0.1974, 0.1974], '--', color='grey', linewidth=1.5)


	## add top axis annotations
	ymin, ymax = ax.get_ylim()
	# ax.text(1.3, ymax-0.1*(ymax-ymin), "%s" % title, fontsize=21)
	# ax.text(1.3, ymax-0.2*(ymax-ymin), r'$\phi$ = %0.2f' % aridity, fontsize=21)

	ax_loc = [7, 30, 90, 180, 365]
	ax_ticks = ['7', '30', '90', '180', '365']
	[ax.text(ax_loc[i], ymax * 1.02, ax_ticks[i], ha='center', fontsize=16) for i in xrange(5)]

	## add yaxis labels and change font size
	ax.set_xlabel(r"$\tau$ (day)", fontsize=18)
	ax.xaxis.set_label_coords(1.05, -0.01)
	ax.set_ylabel("MSC", fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)

	## add legend
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width*0.8, box.height])

	ax.legend(bbox_to_anchor=(1.35, 0.5), loc='center right', fontsize=18)
	# plt.grid(True)

	return

def CoherenceInterVarPlot(ax, coh, index, aridity, fs, freq, title):

	col = ['Crimson', 'Orange', 'DeepSkyblue', 'SpringGreen', 'Purple']
	drivername = ['Tair','Rn','Wind','Humidity','VPD']
	[ax.semilogx(fs/freq, coh[j, :], color=col[i], linewidth=2.5, alpha=0.6, label=drivername[i]) for j, i in enumerate(index)]

	## xlim, ylim, annotation
	ax.set_xlim([1, 2000]) # so that no space before and after the time series

	## add vertical lines
	ax.plot([365, 365], [0, 1], 'r--', linewidth=3.5)
	ax.plot([90, 90], [0, 1], 'b--', linewidth=3.5)
	ax.plot([30, 30], [0, 1], 'k--', linewidth=3.5)
	ax.plot([7, 7], [0, 1], 'g--', linewidth=3.5)
	ax.plot([1, 2000], [0.1974, 0.1974], '--', color='grey', linewidth=1.5)


	## add top axis annotations
	ymin, ymax = ax.get_ylim()
	ax.text(1.3, ymax-0.1*(ymax-ymin), "%s" % title, fontsize=21)

	ax_loc = [7, 30, 90, 365]
	ax_ticks = ['W', 'M', 'S', 'A']
	[ax.text(ax_loc[i], ymax * 1.02, ax_ticks[i], ha='center', fontsize=16) for i in xrange(4)]

	## add yaxis labels and change font size
	ax.set_xlabel(r"$\tau$ (day)", fontsize=18)
	ax.xaxis.set_label_coords(1.05, -0.01)
	ax.set_ylabel("MSC", fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)

	## add legend
	# box = ax.get_position()
	# ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
	#
	# ax.legend(bbox_to_anchor=(1.35, 0.5), loc='center right', fontsize=18)
	plt.grid(True)

	return


def CoherenceModelPlot(ax, coh, aridity, fs, freq, title):

	# col = ['Black', 'Blue', 'Purple', 'Orange', 'Crimson', 'Springgreen']
	# names = ['Pan', 'Penpan', 'Penman', 'PT', 'Hamon', 'Turc']
	nvar = len(coh)
	col = ['Blue', 'Purple', 'Orange', 'Crimson', 'Springgreen']
	names = ['Penpan', 'Penman', 'PT', 'Hamon', 'Turc']
	[ax.semilogx(fs/freq, coh[i, :], marker='o', markersize=4, markeredgecolor=col[i], alpha=0.3, color=col[i], linewidth=2, label=names[i]) for i in xrange(0, nvar)]

	## xlim, ylim, annotation
	ax.set_xlim([1, 2000]) # so that no space before and after the time series
	ymin, ymax = ax.get_ylim()
	# ax.text(1.3, ymax-0.1*(ymax-ymin), "%s" % title, fontsize=21)
	# ax.text(1.3, ymax-0.2*(ymax-ymin), r'$\phi$ = %0.2f' % aridity, fontsize=21)

	## add vertical lines
	ax.plot([365, 365], [0, 1], 'r--', linewidth=3)
	ax.plot([180, 180], [0, 1], 'k--', linewidth=3)
	ax.plot([90, 90], [0, 1], 'g--', linewidth=3)
	ax.plot([30, 30], [0, 1], 'b--', linewidth=3)
	ax.plot([7, 7], [0, 1], 'c--', linewidth=3)
	ax.plot([1, 2000], [0.1974, 0.1974], '--', color='grey', linewidth=1.5)


	## add top axis annotations
	ax_loc = [7, 30, 90, 180, 365]
	ax_ticks = ['7', '30', '90', '180', '365']
	[ax.text(ax_loc[i], ymax * 1.02, ax_ticks[i], ha='center', fontsize=16) for i in xrange(5)]

	## add yaxis labels and change font size
	ax.set_xlabel(r"$\tau$ (day)", fontsize=18)
	ax.xaxis.set_label_coords(1.05, -0.01)
	ax.set_ylabel("MSC", fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)

	## add legend
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width*0.8, box.height])

	ax.legend(bbox_to_anchor=(1.35, 0.5), loc='center right', fontsize=18)
	# plt.grid(True)

	return

def CoherenceModelPointPlot(ax, coh):
	"Plot one variable at a time, for different aridity ranges"
	nvar = 5
	col = ['Blue', 'Purple', 'Orange', 'Crimson', 'Springgreen']
	names = ['Penpan', 'Penman', 'PT', 'Hamon', 'Turc']

	[ax.plot(mean(coh[:, :, i], axis=0), '--s', color=col[i], ms=10, mew=0, linewidth=1.5, alpha=0.6, label=names[i]) for i in xrange(0, nvar)]

	## xlim, ylim, annotation
	ax.set_xlim([-0.5, 5.5]) # so that no space before and after the time series
	ax.set_ylim([0, 1.1])
	## add confidence level
	ax.plot([-1, 6], [0.1974, 0.1974], '--', color='grey', linewidth=1)
	## add confidence level
	ax.plot([-1, 6], [1.0, 1.0], '--', color='black', linewidth=1)

	## add annotations
	# ymin, ymax = ax.get_ylim()

	## add yaxis labels and change font size
	ax.set_xlabel("Period", fontsize=16)
	ax.xaxis.set_label_coords(1.03, -0.02)
	ax.set_ylabel("MSC", fontsize=18)
	ax.set_xticks([0, 1, 2, 3, 4, 5])
	ax.set_xticklabels([7, 30, 90, 120, 180, 365])
	ax.tick_params(axis='x', labelsize=16)
	ax.tick_params(axis='y', labelsize=16)

	## add legend
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
	ax.legend(bbox_to_anchor=(1.35, 0.5), loc='center right', fontsize=16)

	return



def CoherenceEnsemblePlot(ax, coh, fs, freq, title):

	# col = ['Black', 'Blue', 'Purple', 'Orange', 'Crimson', 'Springgreen']
	# names = ['Pan', 'Penpan', 'Penman', 'PT', 'Hamon', 'Turc']
	nvar = len(coh)
	colors = r_[linspace(0, 1, 10), linspace(0, 1, 10)]
	mymap = plt.get_cmap("jet")
	my_colors = mymap(colors)
	[ax.semilogx(fs/freq, coh[i, :], marker='o', markersize=4, markeredgecolor=my_colors[i], alpha=0.3, color=my_colors[i], linewidth=2) for i in xrange(0, nvar)]

	## xlim, ylim, annotation
	ax.set_xlim([1, 2000]) # so that no space before and after the time series
	ymin, ymax = ax.get_ylim()
	# ax.text(1.3, ymax-0.1*(ymax-ymin), "%s" % title, fontsize=21)
	# ax.text(1.3, ymax-0.2*(ymax-ymin), r'$\phi$ = %0.2f' % aridity, fontsize=21)

	## add vertical lines
	ax.plot([365, 365], [0, 1], 'r--', linewidth=3)
	ax.plot([180, 180], [0, 1], 'k--', linewidth=3)
	ax.plot([90, 90], [0, 1], 'g--', linewidth=3)
	ax.plot([30, 30], [0, 1], 'b--', linewidth=3)
	ax.plot([7, 7], [0, 1], 'c--', linewidth=3)
	ax.plot([1, 2000], [0.1974, 0.1974], '--', color='grey', linewidth=1.5)


	## add top axis annotations
	ax_loc = [7, 30, 90, 180, 365]
	ax_ticks = ['7', '30', '90', '180', '365']
	[ax.text(ax_loc[i], ymax * 1.02, ax_ticks[i], ha='center', fontsize=16) for i in xrange(5)]

	## add yaxis labels and change font size
	ax.set_xlabel(r"$\tau$ (day)", fontsize=18)
	ax.xaxis.set_label_coords(1.05, -0.01)
	ax.set_ylabel("MSC", fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)

	## add legend
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width*0.8, box.height])

	# ax.legend(bbox_to_anchor=(1.35, 0.5), loc='center right', fontsize=18)
	plt.grid(True)

	return


def CoherenceWindowPlot(ax, coh, varname, fs, freq):

	nvar = len(coh)
	labels = ['7d', '15d', '30d']
	alps = [0.6, 0.5, 0.4]
	colors = ['green', 'blue', 'red']

	# labels = ['7d', '30d', '90d']
	# colors = ['c', 'blue', 'green']

	[ax.semilogx(fs/freq, coh[i, :], marker='o', markersize=4, markeredgecolor=colors[i], alpha=alps[i], color=colors[i], label=labels[i], linewidth=2) for i in xrange(0, nvar)]

	## xlim, ylim, annotation
	# ax.set_ylim([1, 5*10**6])
	ax.set_ylim([0, 1])
	ax.set_xlim([1, 2000]) # so that no space before and after the time series
	# ax.set_xlim([1, 180]) # so that no space before and after the time series
	ymin, ymax = ax.get_ylim()
	# ax.text(1.3, ymax-0.1*(ymax-ymin), "%s" % title, fontsize=21)
	# ax.text(1.3, ymax-0.2*(ymax-ymin), r'$\phi$ = %0.2f' % aridity, fontsize=21)

	## add vertical lines
	# ax.plot([365, 365], [0, 1], 'c--', linewidth=3)
	ax.plot([90, 90], [0, 1], 'k--', linewidth=3)
	ax.plot([30, 30], [0, 1], 'r--', linewidth=3)
	ax.plot([15, 15], [0, 1], 'b--', linewidth=3)
	ax.plot([7, 7], [0, 1], 'g--', linewidth=3)


	## add top axis annotations
	xmin, xmax = ax.get_xlim()
	ax_loc = [7, 15, 30, 90, 365]
	ax_ticks = ['7', '15', '30', '90', '365']
	[ax.text(ax_loc[i], ymax * 1.02, ax_ticks[i], ha='center', fontsize=14) for i in xrange(5)]
	# ax_loc = [0.2, 0.9]
	# ax_ticks = ['0.2', '0.9']
	# [ax.text(1000, ax_loc[i], ax_ticks[i], va='center', fontsize=15) for i in xrange(2)]
	ax.text(1400, 0.08, "%s" %varname, ha='right', fontsize=18)

	## add yaxis labels and change font size
	ax.set_xlabel(r"$\tau$ (day)", fontsize=14)
	ax.xaxis.set_label_coords(1.05, -0.01)
	ax.set_ylabel("MSC Losses", fontsize=18)
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)

	## add legend
	# box = ax.get_position()
	# ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
	#
	# ax.legend(bbox_to_anchor=(1.35, 0.5), loc='center right', fontsize=18)
	# plt.grid(True)

	return

def CoherenceBasinPlot_3var(ax, coh, aridity, fs, freq, title):

	nvar = coh.shape[0]
	# col = ['LightCoral', 'Orange', 'MediumSpringGreen', 'Orchid', 'SkyBlue', 'Crimson', 'Purple', 'Aquamarine']
	# drivername = ['Tair','Rs','Wind','Humidity','Cloudiness', 'Tsurf', 'VPD', 'Rain']
	col = ['Orange', 'Crimson', 'SpringGreen']
	drivername = ['Rn','Tair','Humidity']
	# linestyle = ['-o','-s','-*','-^','-','-o','-s','-*']
	[ax.semilogx(fs/freq, coh[i, :], color=col[i], linewidth=2.5, alpha=0.6, label=drivername[i]) for i in xrange(nvar-1, -1, -1)]
	# ax.xaxis.set_visible(False)
	plt.xlim([1, 2000]) # so that no space before and after the time series
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)
	ymin, ymax = ax.get_ylim()
	ax.text(1.5, ymax-0.12*(ymax-ymin), "%s" % title, fontsize=20)
	ax.text(1.5, ymax-0.2*(ymax-ymin), r'$\phi$ = %0.2f' % aridity, fontsize=18)
	ax.plot([365, 365], [0, 1], 'r--', linewidth=3.5)
	ax.plot([90, 90], [0, 1], 'k--', linewidth=3.5)
	ax.plot([30, 30], [0, 1], 'b--', linewidth=3.5)
	plt.grid(True)

	return

def CoherenceAridityPlot(ax, coh, index_DI, fs, freq, label, title):
	"Plot one variable at a time, for different aridity ranges"
	# arid_class = [r'0<$\phi$<1 (humid)', r'1<$\phi$<1.5 (sub-humid)', r'1.5<$\phi$<2.5 (transitional)', r'2.5<$\phi$<5 (semi-arid)', r'$\phi$>5 (arid)']
	arid_class = [r'humid (0<$\phi$<2)', r'sub-humid (2<$\phi$<4)', r'semi-arid (4<$\phi$<8)', r'arid ($\phi$>8)']

	colors = r_[linspace(0, 1, 4), linspace(0, 1, 4)]
	mymap = plt.get_cmap("coolwarm")
	my_colors = mymap(colors)
	[ax.semilogx(fs/freq, mean(coh[index, :], axis=0), color=my_colors[i], linewidth=2.5, alpha=0.6, label=arid_class[i]) for i,index in enumerate(index_DI)]

	## xlim, ylim, annotation
	ax.set_xlim([1, 2000]) # so that no space before and after the time series
	ax.set_ylim([0, 1]) # so that no space before and after the time series

	# ax.set_title(label, fontsize=14, y=1.11)

	## add vertical lines, linewidth from 3 to 1
	ax.plot([365, 365], [0, 1], 'r--', linewidth=1.5)
	ax.plot([180, 180], [0, 1], 'k--', linewidth=1.5)
	ax.plot([90, 90], [0, 1], 'g--', linewidth=1.5)
	ax.plot([30, 30], [0, 1], 'b--', linewidth=1.5)
	ax.plot([7, 7], [0, 1], 'c--', linewidth=1)
	ax.plot([1, 2000], [0.1974, 0.1974], '--', color='grey', linewidth=1.5)

	## add annotations
	ymin, ymax = ax.get_ylim()
	# ax.text(1750, 0.08, "(%s) %s" % (label, title), ha='right', fontsize=14)
	ax.text(1.5, 0.84, "(%s)" % (label), ha='left', fontsize=18)
	ax.text(3.5, 0.84, "%s" % (title), ha='left', fontsize=22)
	# ax.text(2, 0.84, "%s" % title, ha='left', fontsize=24)

	# ax.text(1.2, ymax-0.15*(ymax-ymin), "(%s) %s" % (label, title), fontsize=25)
	ax_loc = [7, 30, 90, 180, 365]
	ax_ticks = ['7', '30', '90', '180', '365']
	[ax.text(ax_loc[i], ymax * 1.02, ax_ticks[i], ha='center', fontsize=12) for i in xrange(5)]  # 12

	## add yaxis labels and change font size
	ax.set_xlabel(r"$\tau$ (day)", fontsize=16)
	ax.xaxis.set_label_coords(1.05, -0.01)
	ax.tick_params(axis='x', labelsize=14)

	ax.set_ylabel("MSC", fontsize=18)
	# ax.set_ylabel("Partial Coherence", fontsize=24)
	ax.tick_params(axis='y', labelsize=14)

	# adjust box and grid
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width, box.height])
	# ax.grid(True)

	return

def CoherenceWindow2Plot(ax, coh, fs, freq, label):

	nvar = len(coh)
	labels = ['90d', '30d', '7d']
	alps = [0.5, 0.6, 0.5]
	colors = ['blue', 'green', 'red']

	[ax.semilogx(fs/freq, coh[i, :], marker='o', markersize=4, markeredgecolor=colors[i], alpha=alps[i], color=colors[i], label=labels[i], linewidth=2) for i in xrange(0, nvar)]

	## xlim, ylim, annotation
	ax.set_xlim([1, 2000]) # so that no space before and after the time series
	ax.set_ylim([0, 1]) # so that no space before and after the time series
	ax.set_title(label, fontsize=14, y=1.11)


	## add vertical lines
	ax.plot([365, 365], [0, 1], 'r--', linewidth=3)
	ax.plot([180, 180], [0, 1], 'k--', linewidth=3)
	ax.plot([90, 90], [0, 1], 'g--', linewidth=3)
	ax.plot([30, 30], [0, 1], 'b--', linewidth=3)
	ax.plot([7, 7], [0, 1], 'c--', linewidth=3)
	ax.plot([1, 2000], [0.1974, 0.1974], '--', color='grey', linewidth=1.5)


	## add top axis annotations
	ymin, ymax = ax.get_ylim()
	# ax.text(1.2, ymax-0.15*(ymax-ymin), "(%s) %s" % (label, title), fontsize=25)
	ax_loc = [7, 30, 90, 180, 365]
	ax_ticks = ['7', '30', '90', '180', '365']
	[ax.text(ax_loc[i], ymax * 1.02, ax_ticks[i], ha='center', fontsize=14) for i in xrange(5)]

	## add yaxis labels and change font size
	ax.set_xlabel(r"$\tau$ (day)", fontsize=16)
	ax.xaxis.set_label_coords(1.05, -0.01)
	ax.set_ylabel("MSC", fontsize=18)
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)

	## add legend
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width, box.height])
	ax.grid(True)

	return



def CoherenceWindowPointPlot(ax, coh, label):

	"Plot one variable at a time, for different aridity ranges"
	nvar = coh.shape[1]
	labels = ['90d', '30d', '7d']
	alps = [0.5, 0.6, 0.5]
	colors = ['blue', 'green', 'red']

	[ax.plot(coh[:, i], '--o', color=colors[i], ms=10, mew=0, linewidth=1.5, alpha=alps[i], label=labels[i]) for i in xrange(0, nvar)]

	## xlim, ylim, annotation
	ax.set_xlim([-0.5, 5.5]) # so that no space before and after the time series
	ax.set_ylim([0, 1.1])
	ax.set_title(label, fontsize=16, y=1.02)  #  y=1.11

	## add annotations
	# ymin, ymax = ax.get_ylim()

	## add yaxis labels and change font size
	ax.set_xlabel("Period", fontsize=14)
	ax.xaxis.set_label_coords(1.08, -0.02)
	ax.set_ylabel("MSC Losses", fontsize=18)
	ax.set_xticks([0, 1, 2, 3, 4, 5])
	ax.set_xticklabels([7, 30, 90, 120, 180, 365])
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)

	return

def CoherenceWindowAridityPointPlot(ax, coh, index_DI, label):
	"Plot one variable at a time, for different aridity ranges"
	arid_class = [r'humid (0<$\phi$<2)', r'sub-humid (2<$\phi$<4)', r'semi-arid (4<$\phi$<8)', r'arid ($\phi$>8)']
	colors = r_[linspace(0, 1, 4), linspace(0, 1, 4)]
	mymap = plt.get_cmap("coolwarm")
	my_colors = mymap(colors)
	[ax.plot(mean(coh[:, index], axis=1), '--o', color=my_colors[i], ms=10, mew=0, linewidth=1.5, alpha=0.6, label=arid_class[i]) for i,index in enumerate(index_DI)]

	## xlim, ylim, annotation
	ax.set_xlim([-0.5, 5.5]) # so that no space before and after the time series
	ax.set_ylim([0, 1.1])
	ax.set_title(label, fontsize=16, y=1.02)  #  y=1.11

	## add yaxis labels and change font size
	ax.set_xlabel("Period", fontsize=14)
	ax.xaxis.set_label_coords(1.08, -0.02)
	ax.set_ylabel("MSC Losses", fontsize=18)
	ax.set_xticks([0, 1, 2, 3, 4, 5])
	ax.set_xticklabels([7, 30, 90, 120, 180, 365])
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)

	return

def CoherencePointAridityPlot(ax, coh, index_DI, label, title):
	"Plot one variable at a time, for different aridity ranges"
	arid_class = [r'humid (0<$\phi$<2)', r'sub-humid (2<$\phi$<4)', r'semi-arid (4<$\phi$<8)', r'arid ($\phi$>8)']
	colors = r_[linspace(0, 1, 4), linspace(0, 1, 4)]
	mymap = plt.get_cmap("coolwarm")
	my_colors = mymap(colors)
	[ax.plot(mean(coh[index, :], axis=0), '--o', color=my_colors[i], ms=12, mew=0, linewidth=1.5, alpha=0.6, label=arid_class[i]) for i,index in enumerate(index_DI)]
	# xloc = array([7, 30, 90, 180, 365])
	# [ax.plot(xloc-1, mean(coh[index, :], axis=0), '--o', color=my_colors[i], ms=12, mew=0, linewidth=1.5, alpha=0.6, label=arid_class[i]) for i,index in enumerate(index_DI)]

	## xlim, ylim, annotation
	ax.set_xlim([-1, 5]) # so that no space before and after the time series
	ax.set_ylim([0, 1.1])
	## add confidence level
	ax.plot([-1, 6], [0.1974, 0.1974], '--', color='grey', linewidth=1)

	## add annotations
	ymin, ymax = ax.get_ylim()
	ax.text(-0.7, ymax-0.18*(ymax-ymin), "(%s) %s" % (label, title), fontsize=18)

	## add yaxis labels and change font size
	ax.set_xlabel("Period", fontsize=14)
	ax.xaxis.set_label_coords(1.02, -0.02)
	ax.set_ylabel("MSC", fontsize=18)
	ax.set_xticks([0, 1, 2, 3, 4])
	ax.set_xticklabels([7, 30, 90, 180, 365])
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)

	return

def CoherencePointSingleAridityPlot(ax, coh, index_DI, index_model, label, title):
	"Plot one variable at a time, for different aridity ranges"
	models = ['Epan', 'Penpan', 'Penman', 'PT', 'Hamon', 'Turc']
	col = ['Black', 'Blue', 'Purple', 'Orange', 'Crimson', 'Springgreen']
	ms = [11, 7, 7, 7, 7, 7]  # 8.5
	linestyle = ['o', 's', 's', 's', 's', 's']
	alp = [0.55, 0.5, 0.5, 0.5, 0.5, 0.5]
	monmsc = mean(coh[index_DI, :, :], axis=0)
	[ax.plot(monmsc[:, i], '--', color=col[i], marker=linestyle[i], ms=ms[i], mew=0, linewidth=1.5, alpha=alp[i], label=models[i]) for i in index_model]

	## xlim, ylim, annotation
	ax.set_xlim([-0.5, 5.5]) # so that no space before and after the time series
	ax.set_ylim([0, 1.1])
	## add confidence level
	ax.plot([-1, 6], [0.1974, 0.1974], '--', color='grey', linewidth=1)

	## add annotations
	ymin, ymax = ax.get_ylim()
	ax.text(5.3, ymin+0.19*(ymax-ymin), "(%s) %s" % (label, title), ha='right', fontsize=15)

	## add yaxis labels and change font size
	ax.set_xlabel("Period", fontsize=14)
	ax.xaxis.set_label_coords(1.1, -0.03)
	ax.set_ylabel("MSC", fontsize=18)
	ax.set_xticks([0, 1, 2, 3, 4, 5])
	ax.set_xticklabels([7, 30, 90, 120, 180, 365])
	ax.tick_params(axis='x', labelsize=13)
	ax.tick_params(axis='y', labelsize=13)

	return

def ModelErrorBarPlot(ax, stat):
	"Plot one variable at a time, for different aridity ranges"
	models = ['Penpan', 'Penman', 'PT', 'Hamon', 'Turc']
	col = ['Blue', 'Purple', 'Orange', 'Crimson', 'Springgreen']
	ms = [8.5, 8.5, 8.5, 8.5, 8.5]
	linestyle = ['s', 's', 's', 's', 's']
	alp = [0.5, 0.5, 0.5, 0.5, 0.5]

	[ax.errorbar(np.arange(4), stat[:, 0, i],  stat[:, 1, i], color=col[i], linewidth=2, capsize=10, elinewidth=5) for i, model in enumerate(models)] #  alpha=alp[i], label=model, marker=linestyle[i], ms=ms[i], mew=0,

	## xlim, ylim, annotation
	ax.set_xlim([-1, 4]) # so that no space before and after the time series

	# ## add annotations
	# ymin, ymax = ax.get_ylim()
	# ax.text(-0.7, ymax-0.18*(ymax-ymin), "(%s) %s" % (label, title), fontsize=18)
	#
	# ## add yaxis labels and change font size
	# ax.xaxis.set_label_coords(1.02, -0.02)
	ax.set_ylabel("Residue", fontsize=18)
	ax.set_xticks([0, 1, 2, 3])
	ax.set_xticklabels(['Daily', 'Weekly', 'Monthly', 'Annual'])
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)

	return

def CoherenceSingleTimePointModelPlot(ax, coh, index_DI, index_model, label, title):
	"Plot one variable at a time, for different aridity ranges"
	models = ['Pan', 'Penpan', 'Penman', 'PT', 'Hamon', 'Turc']
	col = ['Black', 'Blue', 'Purple', 'Orange', 'Crimson', 'Springgreen']
	ms = [12, 8.5, 8.5, 8.5, 8.5, 8.5]
	linestyle = ['o', 's', 's', 's', 's', 's']
	alp = [0.4, 0.5, 0.5, 0.5, 0.5, 0.5]

	monmsc = array([mean(coh[index, :], axis=0) for index in index_DI])
	[ax.plot(monmsc[:, i], '--', color=col[i], marker=linestyle[i], ms=ms[i], mew=0, linewidth=1.5, alpha=alp[i], label=models[i]) for i in index_model]

	## xlim, ylim, annotation
	ax.set_xlim([-1, 4]) # so that no space before and after the time series
	ax.set_ylim([0, 1.1])
	## add confidence level
	ax.plot([-1, 5], [0.1974, 0.1974], '--', color='grey', linewidth=1)

	## add annotations
	ymin, ymax = ax.get_ylim()
	ax.text(-0.7, ymax-0.15*(ymax-ymin), "(%s) %s" % (label, title), fontsize=14)

	## add yaxis labels and change font size
	ax.set_xlabel("Dryness", fontsize=14)
	ax.xaxis.set_label_coords(1.05, -0.02)
	ax.set_ylabel("MSC", fontsize=18)
	ax.set_xticks([0, 1, 2, 3])
	ax.set_xticklabels(['humid', 'sub-humid', 'semi-arid', 'arid'], rotation=15)
	ax.tick_params(axis='x', labelsize=13)
	ax.tick_params(axis='y', labelsize=14)

	return



def CoherenceBasinPlotTogether(ax, coh, fs, freq, basinnames):

	nvar = coh.shape[0]
	# col = ['LightCoral', 'Orange', 'MediumSpringGreen', 'Orchid', 'SkyBlue', 'Crimson', 'Purple', 'Aquamarine', 'red', 'black']
	cm_sub = linspace(0, 1, 10)
	col = [plt.cm.gist_rainbow(x) for x in cm_sub]
	# linestyle = ['-o','-s','-*','-^','-','-o','-s','-*']
	[ax.semilogx(fs/freq, coh[i, :], color=col[i], linewidth=5, alpha=0.2, label=basinnames[i]) for i in xrange(0, nvar)]
	# ax.xaxis.set_visible(False)
	plt.xlim([1, 2000]) # so that no space before and after the time series
	# ax.set_xticks(years[2::20])  # specify where you put the ticks
	# ax.set_xticklabels(years[2::20], fontsize=15)  # adjust the ticks font size
	# ax.set_ylabel("Z-score", fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)
	# ax.set_title("%s" % title, fontsize=18, y=1.02)
	# ymin, ymax = ax.get_ylim()
	# ax.text(1.5, ymax-0.15*(ymax-ymin), "%s" % title, fontsize=20)
	ax.plot([365, 365], [0, 1], 'r--', linewidth=3.5)
	ax.plot([90, 90], [0, 1], 'k--', linewidth=3.5)
	ax.plot([30, 30], [0, 1], 'b--', linewidth=3.5)
	plt.grid(True)

	return

def Coherence_obs_mod_Plot(ax, coh, fs, freq, title):

	# linestyle = ['-o','-s','-*','-^','-','-o','-s','-*']
	ax.semilogx(fs/freq, coh, color='grey', linewidth=3)
	# ax.xaxis.set_visible(False)
	plt.xlim([1, 2000]) # so that no space before and after the time series
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)
	ymin, ymax = ax.get_ylim()
	ax.text(1.5, ymax-0.15*(ymax-ymin), "%s" % title, fontsize=20)
	ax.plot([365, 365], [0, 1], 'r--', linewidth=3.5)
	ax.plot([90, 90], [0, 1], 'k--', linewidth=3.5)
	ax.plot([30, 30], [0, 1], 'b--', linewidth=3.5)
	plt.grid(True)

	return

def Mapshow(data, lons, lats, size, alpha, min, max, cmp, title, legend, unit, figdir, filename):

	m = Basemap(projection='merc', llcrnrlon=71, llcrnrlat=3, urcrnrlon=139, urcrnrlat=55, lat_ts=20)
	# draw boundaries using shapefile with countries boundaries
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou1_4m/bou1_4l', 'CHN_adm1', linewidth=2)
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou2_4m/bou2_4l', 'CHN_adm2', linewidth=1)
	m.drawparallels(np.arange(10, 60, 20), labels=[1, 0, 0, 0], fontsize=19)  # only left ytick
	m.drawmeridians(np.arange(90, 140, 30), labels=[0, 0, 0, 1], fontsize=19)  # only bottom xtick

	# draw data
	if type(data) is str:
		im = m.scatter(lons, lats, size, marker='o', c=data, latlon=True, label=legend)
	else:
		im = m.scatter(lons, lats, size, marker='o', c=data, edgecolors='none', alpha=alpha, vmin=min, vmax=max, latlon=True, cmap=cmp, label=legend)
		cb = m.colorbar(im, pad='3%')
		cb.ax.tick_params(labelsize=20)
	plt.title(title, fontsize=16) #, y=1.02)
	# plt.legend(loc=3, prop={'size': 18})
	plt.xlabel(unit, fontsize=21, labelpad=19)
	# savefig('%s%s' % (figdir, filename))
	# plt.show()
	# plt.clf()
	return

def Mapshow_Aridity(lons, lats, size, alpha, min, max, color, legend):

	m = Basemap(projection='merc', llcrnrlon=71, llcrnrlat=3, urcrnrlon=139, urcrnrlat=55, lat_ts=20)
	# draw boundaries using shapefile with countries boundaries
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou1_4m/bou1_4l', 'CHN_adm1', linewidth=0.9)
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou2_4m/bou2_4l', 'CHN_adm2', linewidth=0.5)
	m.drawparallels(np.arange(10, 60, 20), labels=[1, 0, 0, 0], fontsize=14, linewidth=0.0)  # only left ytick
	m.drawmeridians(np.arange(90, 140, 30), labels=[0, 0, 0, 1], fontsize=14, linewidth=0.0)  # only bottom xtick

	m.scatter(lons, lats, size, marker='o', c=color, edgecolors='none', alpha=alpha, vmin=min, vmax=max, latlon=True, label=legend, rasterized=True)
	plt.legend(loc=3, prop={'size': 14})
	plt.hold(True)
	return

def Mapshow_Aridity_subplot(ax, lons, lats, size, alpha, min, max, color, legend):

	ax.m = Basemap(projection='merc', llcrnrlon=71, llcrnrlat=3, urcrnrlon=139, urcrnrlat=55, lat_ts=20)
	# draw boundaries using shapefile with countries boundaries
	ax.m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou1_4m/bou1_4l', 'CHN_adm1', linewidth=0.9)
	ax.m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou2_4m/bou2_4l', 'CHN_adm2', linewidth=0.5)
	ax.m.drawparallels(np.arange(10, 60, 20), labels=[1, 0, 0, 0], fontsize=14, linewidth=0.0)  # only left ytick
	ax.m.drawmeridians(np.arange(90, 140, 30), labels=[0, 0, 0, 1], fontsize=14, linewidth=0.0)  # only bottom xtick

	ax.m.scatter(lons, lats, size, marker='o', c=color, edgecolors='none', alpha=alpha, vmin=min, vmax=max, latlon=True, label=legend)
	plt.annotate('(a)', xy=(0.08,0.87), xycoords='axes fraction', fontsize=16)
	plt.legend(loc=3, prop={'size': 10})
	plt.hold(True)
	return

def Mapshow_RGB(data, lons, lats, size, alpha, marker, title, legend, unit):

	m = Basemap(projection='merc', llcrnrlon=71, llcrnrlat=3, urcrnrlon=139, urcrnrlat=55, lat_ts=20)
	# draw boundaries using shapefile with countries boundaries
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou1_4m/bou1_4l', 'CHN_adm1', linewidth=2)
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou2_4m/bou2_4l', 'CHN_adm2', linewidth=1)
	m.drawparallels(np.arange(10, 60, 20), labels=[1, 0, 0, 0], fontsize=19)  # only left ytick
	m.drawmeridians(np.arange(90, 140, 30), labels=[0, 0, 0, 1], fontsize=19)  # only bottom xtick
	# draw data
	im = m.scatter(lons, lats, size, marker=marker, c=data, edgecolors=data, alpha=alpha, latlon=True, label=legend)
	# cb = m.colorbar(im, pad='3%')
	# plotdata = m.transform_scalar(data, lons, lats, nx, ny)
	# im = m.imshow(plotdata, vmin=min, vmax=max, cmap=cmp)
	plt.title(title, fontsize=16, y=1.02)
	plt.xlabel(unit, labelpad=20)
	# plt.show()
	# plt.clf()

def MapShow_lambert():
	m = Basemap(llcrnrlon=82, llcrnrlat=0, urcrnrlon=140, urcrnrlat=55, projection='lcc', lat_1=20, lat_2=40, lon_0=108)
	# draw boundaries using shapefile with countries boundaries
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou1_4m/bou1_4l', 'CHN_adm1', linewidth=2)
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou2_4m/bou2_4l', 'CHN_adm2', linewidth=1)
	# m.etopo()
	m.drawparallels(np.arange(20, 60, 20), labels=[1, 0, 0, 0])  # only left ytick
	m.drawmeridians(np.arange(90, 140, 30), labels=[0, 0, 0, 1])  # only bottom xtick
	plt.show()


def TimeseriesPlot_vars(input, styr, edyr, title):

	fig = plt.figure(figsize=(14, 8))
	ax = fig.add_subplot(1, 1, 1)
	years = np.arange(styr, edyr+1)

	col = ['Crimson', 'Orange', 'DeepSkyblue', 'SpringGreen', 'Purple', 'Black']
	drivername = ['Tair','Rn','Wind','Humidity','VPD', 'Pan']
	widths = [3.5, 3.5, 3.5, 3.5, 3.5, 4]
	# linestyle = ['-o','-s','-*','-^','-', '-']

	for i in [5, 4, 1, 0, 3, 2]:
		ts = (input[i,:] - mean(input[i, :]))/std(input[i, :])
		ts_rolling = pd.rolling_mean(ts, 5)
		R = corrcoef(input[5, :], input[i,:])[0,1]
		if i < 5:
			ax.plot(years, ts_rolling, color=col[i], label="%s (%.02f)" % (drivername[i], R**2), linewidth=widths[i], alpha=0.6)  #linestyle[i],
		else:
			ax.plot(years, ts_rolling, color=col[i], label="%s" % (drivername[i]), linewidth=widths[i], alpha=0.6)  #linestyle[i],


	plt.xlim([styr, edyr]) # so that no space before and after the time series
	ax.set_xticks(years[4::5])  # specify where you put the ticks
	ax.set_xticklabels(years[4::5], fontsize=15)  # adjust the ticks font size
	ax.set_ylabel("Z-score", fontsize=24)
	ax.tick_params(axis='y', labelsize=18)
	plt.legend(loc=0)
	ymin, ymax = ax.get_ylim()
	ax.text(1965, ymax-0.2*(ymax-ymin), "%s" % title.upper(), fontsize=24)
	plt.show()

	return