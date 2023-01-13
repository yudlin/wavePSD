def wavePSD(X, fs, flo, fup, fres, sFig=True):
    '''
    This program is used to evaluate the power spectrum of X by complex
    Morlet wavelet.
    Designed by Yue-Der Lin, Feng Chia University, Taiwan.

    Inputs:
           X    The signal to be analyzed.
          fs    Sampling frequency of X.
          flo    Lower bound for spectrum analysis.
         fup    Upper bound for spectrum analysis.
        fres    Frequency resolution in the spectrum.
        sFig    True = Show the figures, False = Do not show the figures
               (default = True).

    Outputs:
          Ta    Times for spectrogram.
       Freqs    Frequencies for spectrogram and spectrum.
         stX    Spectrogram of X.
        psdX    Averaged spectrum of X.

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
    '''
    # Import modules:
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    import pywt

    # Initialization:
    matplotlib.rcParams['font.family'] = 'Times'
    #   For complex Morlet wavelet:
    #     fs = fc/delta, fc = 1 for 'cmor1-1'.
    delta = 1/fs
    Fa = np.arange(flo, fup+fres,  fres)  # Frequencies for analysis:
    S = np.divide(np.array([fs]), Fa)
    L = np.size(X)
    Ta = np.arange(0, L/fs, 1/fs)

    #     CoefX, StX: Arrays of (len(S), L).
    CoefX, Freqs = pywt.cwt(X, S, 'cmor1-1', delta)
    StX = np.abs((np.multiply(CoefX, np.matrix.conjugate(CoefX)))/L)
    #     psdX: Array of (lens(s), ).
    psdX = np.mean(StX, axis=1)

    if sFig == True:
        # Spectrogram for X:
        plt.figure(1, figsize = (7.2, 5.4))
        plt.pcolormesh(Ta, Freqs, StX, cmap='ocean_r', shading='nearest')
        plt.xlim(0, L/fs)
        plt.ylim(flo, fup)  # Show spectrogram in frequency range [flo, fup].
        plt.colorbar()
        plt.xlabel('Time (sec)', fontsize = 13)
        plt.ylabel('Frequency (Hz)', fontsize = 13)
        plt.title('Spectrogram', fontsize = 15)
        plt.savefig('Fig_Spectrogram.jpg', dpi=1200, transparent=True)

        # Averaged Spectrum for X:
        plt.figure(2, figsize = (7.2, 5.4))
        plt.plot(Freqs, psdX)
        plt.xlim(flo, fup)  # Show spectrum in frequency range [flo, fup].
        plt.xlabel('Frequency (Hz)', fontsize = 13)
        plt.ylabel('Power', fontsize = 13)
        plt.title('Averaged Spectrum', fontsize = 15)
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
        plt.savefig('Fig_Averaged_Spectrum.jpg', dpi=1200, transparent=True)

        # Combination of Averaged Spectrum and SSQW Spectrogram for X:
        fig3 = plt.figure(3, figsize=(8, 4.2))
        plt.subplot(1,2,2)
        plt.pcolormesh(Ta, Freqs, StX, cmap = 'ocean_r', shading = 'nearest')
        plt.colorbar()
        plt.xlabel('Time (sec)', fontsize = 13)
        # plt.ylabel(‘Frequency (Hz)', fontsize = 13)
        plt.title('Spectrogram', fontsize = 15)
        plt.xlim(0, L/fs); plt.ylim(flo, fup)
        # fig3.subplots_adjust(wspace = 0.4)
        plt.subplot(1,2,1)
        plt.plot(psdX,Freqs)
        plt.xlabel('Power', fontsize = 13)
        plt.ylabel('Frequency (Hz)', fontsize = 13)
        plt.title('Averaged Spectrum', fontsize = 15)
        plt.ylim(flo, fup)
        plt.ticklabel_format(axis='x', style='sci', scilimits=(-3,3))
        # fig3.subplots_adjust(wspace = 0.4)
        plt.savefig('Fig_Wavelet_Spectrum_Spectrogram.jpg', dpi=1200, transparent=True)
    else:
        pass

    return Ta, Freqs, StX, psdX