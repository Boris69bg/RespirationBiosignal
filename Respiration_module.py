# ---------------

# This module provides methods to process Respiration signals, following the approach described by:
# mConverse: Inferring Conversation Episodes fromRespiratory Measurements Collected in the Field.

#----------------

import numpy as np

# local
from biosppy import utils
import matplotlib.pyplot as plt
import h5py
from biosppy.signals import resp



def findQuartiles (signal=None):
	"""
	# auxiliar function for the 1st processing step (outlier removal), by using quartiles
	# Since quartiles are less sensitive to spikes that may appear in respiration measurements, those can be used for outlier detection. 
	# Finding the upper quartile (UQ), the lower quartile (LQ), and the interquartile range (IQR) = UQ-LQ. 
    
    Parameters
    ----------
    signal : array
        Raw Respiration signal.

    Returns
    -------
    UQ :
    	Upper Quartile
	LQ :
		Lower Quartile
	IQR :
		InterQuartile Range

	"""

	# check inputs
	if signal is None:
		raise TypeError("Please specify an input signal.")

	# ensure numpy
	signal = np.array(signal)

	#sort the signal
	sortedS=np.sort(signal)

	# find the upper quartile (UQ=Q4) and lower quartile (LQ=Q1)
	Q1=sortedS[len(signal)/4]
	Q4=sortedS[len(signal)*3/4]

	#compute the interquartile range (IQR) = UQ-LQ
	IQR=Q4-Q1
	# output
	args = (Q1,Q4,IQR)
	names = ('LQ','UQ','IQR')

	return utils.ReturnTuple(args, names)


def outlierdetection (signal=None, UQ=None, LQ=None,IQR=None ):	
	"""
	# 1st Step - Outlier Detection
	# outliers' condition: those points that are >= 1.5IQR+UQ or <= LQ-1.5IQR.

	
    Parameters
    ----------
    signal : array
        Raw Respiration signal.
    UQ, LQ, IQR :
    	Upper and Lower Quartiles, and InterQuartile Range

    Returns
    -------
	outliers : list
		outliers indices
	perc_corrupted : float 
		percentage of outliers/corrupted signal

	"""

	# check inputs
	if signal is None:
		raise TypeError("Please specify an input signal.")

	# ensure numpy
	signal = np.array(signal)	

	LQ,UQ,IQR=findQuartiles(signal)

	outliers=[]

	for q in range (0,len(signal)):
		if signal[q]>=(1.5*IQR+UQ) or signal[q]<=(LQ-1.5*IQR):
			outliers.append(q)

	perc_corrupted=round(float(len(outliers))/float(len(signal))*100,3)

	# output 
	args = (outliers,perc_corrupted)
	names = ('outliers','perc_corrupted')

	return utils.ReturnTuple(args, names)

			

def PeaksID (signal=None,UQ=None, LQ=None,IQR=None):
	"""
	# 2nd Step - Identifying Respiration Cycles by finding its peaks
	# Search of Peaks is computed in 1500 samples (1.5) seconds windows
		
    Parameters
    ----------
    signal : array
        Raw Respiration signal.
    UQ, LQ, IQR :
    	Upper and Lower Quartiles, and InterQuartile Range

    Returns
    -------
	indexPeaks : list
		indices of the Peaks
	peaks : list
		Peaks of the input signal
	"""

	# check inputs
	if signal is None:
		raise TypeError("Please specify an input signal.")

	# ensure numpy
	signal = np.array(signal)	

	LQ,UQ,IQR=findQuartiles(signal)

	indexPeaks=[]
	lens=len(signal)
	peaks=[]
	for t in range(0, len(signal),1500):
		portion=signal[t:t+1500]
		localmax=np.abs(portion.max())
		localmin=np.abs(portion.min())

		if localmax>UQ:
			peaks.append(localmax)
			# if last window has less than 1500 samples
			if t+1500>lens:
				for t in range(t,lens):
					if signal[t]==localmax:
						indexPeaks.append(t)
			else:
				for t in range(t,t+1500):
					if signal[t]==localmax:
						indexPeaks.append(t)
						break

	if (indexPeaks[0]==0):
		indexPeaks= indexPeaks[1:]

	# output
	args = (indexPeaks,peaks)
	names = ('indexPeaks','peaks')

	return utils.ReturnTuple(args, names)



def IRR (signal=None, peaks=None, sampling_rate=1000.,physiological_lim=0.45):
	"""
	# 3rd Step - Selecting only the physiologically meaningful Respiration Cycles and respective Peaks
				
    Parameters
    ----------
    signal : array
        Raw Respiration signal.
    UQ, LQ, IQR :
    	Upper and Lower Quartiles, and InterQuartile Range
	physiological_lim : float, optional
		Only the IRR values below this physiological threshold are accepted.
		Default is 0.45 
		( no biosppy estava 0.35 - confirmar qual o melhor valor para default, nos dados deste dataset 0.45 pareceu razoavel)

    Returns
    -------
	peaks_corrected : list
		indices of the Peaks after physiological-based selection
	rate : list
		instantaneous respiration rate
	"""

	# Original indices of the peaks
	peaks=PeaksID(signal=data5)['indexPeaks']

	# Instantaneous respiration rate
	rate= sampling_rate * (1. / np.diff(peaks)) 

	# Indices of peaks w/ physiological meaning
	peaks_corrected = np.asarray(peaks)					

	# accepting only the IRR within physiological limits
	indx = np.nonzero(rate <= physiological_lim) 
	rate = rate[indx]
	if len(indx[0]>0):
		newindx=np.append(indx[0],indx[0][-1]+1)
		peaks_corrected = peaks_corrected[newindx]
	else:
		peaks_corrected = peaks_corrected[indx]

	# get time vectors
	length = len(signal)
	T = (length - 1) / sampling_rate
	ts = np.linspace(0, T, length, endpoint=False)

	# output
	args = (peaks_corrected, rate)
	names = ('peaks_corrected', 'rate')

	return utils.ReturnTuple(args, names)



# Example:

string="Participant"
hdf5=".hdf5"
filename=string+str(1)+hdf5
f = h5py.File(filename, 'r')

data5=list(f['Video'+str(1)]['1']['Respiration'])
print(outlierdetection(signal=data5)['perc_corrupted'])
print(len(outlierdetection(signal=data5)['outliers']))
'''
ind_Peaks=PeaksID(signal=data5)['indexPeaks']
ind_Peaks_corrected=IRR(signal=data5)['peaks_corrected']

plt.plot(np.array(data5))
for u in range(0, len(ind_Peaks)):
	plt.axvline(ind_Peaks[u], c='r')
plt.show()

plt.plot(np.array(data5))
for u in range(0, len(ind_Peaks_corrected)):	
	plt.axvline(ind_Peaks_corrected[u], c='g')
plt.show()

# with Biosppy:
ind_Peaks=resp.resp(signal=data5,show=False)['zeros']

plt.plot(np.array(data5))
for u in range(0, len(ind_Peaks)):	
	plt.axvline(ind_Peaks[u], c='g')
plt.show()
'''
