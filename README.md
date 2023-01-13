# wavePSD
This program is used to evaluate the power spectrum of X by complex     Morlet wavelet.

Resource:
     
     The PyWavelets package can be found in the following link
     
         https://pywavelets.readthedocs.io/en/latest/install.html

Usage Example:

import numpy as np

sample_rate = 100

seconds = 10

num_samples = sample_rate*seconds

time_vect = np.linspace(0, seconds, num_samples)

x = 0.5*np.cos(2*np.pi*3*time_vect+0.1)

from wavePSD import wavePSD

Ta, Freqs, StX, psdX = wavePSD(x, 100, 1, 10, 0.1, True)
