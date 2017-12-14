import numpy as np
import matplotlib.pyplot as plt

def plomp(f1,f2):
    """ Plomp & Levelt consonance formula for two sinusoids of the same loudness

    Parameters
    ----------
    f1,f2 : frequencies of the two sinusoids


    Returns
    -------
    A consonance value according to the Plomp & Levelt formula
    """
    fmin=min(f1,f2)
    fmax=max(f1,f2)
    s=0.24/(0.021*fmin+19.)
    return (np.exp(-3.5*s*(fmax-fmin))-np.exp(-5.75*s*(fmax-fmin)))


def plomp_spectrum(spectrum):
    """ Plomp & Levelt consonance formula for p spectrums of q harmonics each

    Parameters
    ----------
    spectrum : a p x q numpy array, the (i,j) value of which is the frequency
    of the j-th harmonic of the i-th sound.

    Returns
    -------
    A consonance value for the whole sum of sinusoids
    """
    theShape=spectrum.shape
    nSpectr=theShape[0]
    nFreq=theShape[1]

    c=0.0
    for i in range(nSpectr):
        for j in range(i+1,nSpectr):
            for k in range(nFreq):
                for l in range(nFreq):
                    c=c+plomp(spectrum[i][k],spectrum[j][l])
    return c


########################
###### MAIN

## We will build a consonance table

chordNames=['$C_M$','$C\sharp_M$','$D_M$','$Eb_M$','$E_M$','$F_M$',
            '$F\sharp_M$','$G_M$','$G\sharp_M$','$A_M$','$Bb_M$','$B_M$',
            '$C_m$','$C\sharp_m$','$D_m$','$Eb_m$','$E_m$','$F_m$',
            '$F\sharp_m$','$G_m$','$G\sharp_m$','$A_m$','$Bb_m$','$B_m$']

## Reading the file containing the 96_tunings
## The file contains the number of tunings
## then for each one of them, a line for its name, and the tuning in cents
## for the twelve notes on the following line

tuning_names=[]

with open("./96_tunings.txt","r") as f:
    num_tunings = int(f.readline())
    tunings = np.zeros((num_tunings,12))
    for i in range(num_tunings):
        tuning_names.append(f.readline())
        tunings[i,:] = [float(x) for x in f.readline().rstrip().split(",")]

## Consonance table calculations

consonance_table = np.zeros((tunings.shape[0],24))

for tuning_idx,tuning in enumerate(tunings):
    ## For each tuning, we calculate the corresponding scale tuning_frequencies
    ## assuming a given base frequency

    tuning = list(tuning)
    extended_tuning = tuning + [1200+x for x in tuning]
    base_frequency = 260.
    tuning_frequencies = [base_frequency*np.power(2.0,x/1200.0) for x in extended_tuning]

    ## For reference, we calculate the consonance values of a pure major chord
    ## and a pure minor chord at the given base frequency

    puremajor_spectrum = np.zeros((3,6))
    puremajor_spectrum[:,0] = base_frequency*np.array([1.0,1.25,1.5])
    for i in range(1,6):
        puremajor_spectrum[:,i] = float(i+1)*puremajor_spectrum[:,0]
    puremajor_c = plomp_spectrum(puremajor_spectrum)

    pureminor_spectrum = np.zeros((3,6))
    pureminor_spectrum[:,0] = base_frequency*np.array([1.0,1.2,1.5])
    for i in range(1,6):
        pureminor_spectrum[:,i] = float(i+1)*pureminor_spectrum[:,0]
    pureminor_c = plomp_spectrum(pureminor_spectrum)


    ## For each degree of the scale, we calculate the consonance values of a
    ## major chord and a minor chord which uses the given tuning
    ## Notice that we rescale the frequency ratios to the base base_frequency
    ## since the Plomp & Levelt values decrease with increasing frequency
    ## The consonance table contains the consonance deviations from the pure
    ## major and minor chords

    for index in range(12):
        spectrum = np.zeros((3,6))
        spectrum[:,0] = np.array([tuning_frequencies[index],
                                  tuning_frequencies[index+4],
                                  tuning_frequencies[index+7]])
        for i in range(1,6):
            spectrum[:,i] = float(i+1)*spectrum[:,0]
        spectrum = 260.*spectrum/spectrum[0,0]
        spectrum_c = plomp_spectrum(spectrum)
        consonance_table[tuning_idx,index] = spectrum_c-puremajor_c

        spectrum = np.zeros((3,6))
        spectrum[:,0] = np.array([tuning_frequencies[index],
                                  tuning_frequencies[index+3],
                                  tuning_frequencies[index+7]])
        for i in range(1,6):
            spectrum[:,i] = float(i+1)*spectrum[:,0]
        spectrum = 260.*spectrum/spectrum[0,0]
        spectrum_c = plomp_spectrum(spectrum)
        consonance_table[tuning_idx,index+12] = spectrum_c-pureminor_c

## We display the consonance table in a heatmap
plt.imshow(consonance_table,cmap="magma")
plt.colorbar()
plt.yticks(range(num_tunings),tuning_names,fontsize=6)
plt.xticks(range(24),chordNames,fontsize=6)
plt.show()

## We can also display the deviations relative to the equal temperament
## which is the first tuning in the tuning array
plt.imshow(consonance_table-np.tile(consonance_table[0,:],(num_tunings,1)),cmap="RdYlGn_r")
plt.clim([-0.5,0.5])
plt.colorbar()
plt.yticks(range(num_tunings),tuning_names,fontsize=6)
plt.xticks(range(24),chordNames,fontsize=6)
plt.show()
