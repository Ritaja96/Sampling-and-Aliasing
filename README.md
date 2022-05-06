# Setup
If no version is specified any should work. Here, we used
* MATLAB 2021b
* DSPFirst Toolbox
# How to install and set up DSPFirst toolbox into MATLAB
Install:
Add spfirst_zip files to MATLAB and update MATLAB's path.(link: https://dspfirst.gatech.edu/matlab/spfirst_zip/spfirst_v174.zip)

Steps to follow:
* Download the ZIP file called spfirst_vNNN.zip. that contains everything.
* Unzip spfirst_vNNN.zip and set up like this : \MATLAB\toolbox\spfirst\. The unzip will create a bunch of subdirectories.
**NOTE: use forward slash for directories on Mac or Linux.**
* Then add the above directory path to the MATLAB path. This can be done in MATLAB from the Home Tab ➤ Environment ➤ Set Path.
**NOTE: For versions previous to R2012a, select from menu File ➤ Set Path.**
* Type the command spfirst at the MATLAB comand prompt. It will add the other appropriate subdirectories to the path.
* Final step: under Home Tab ➤ Environment ➤ Set Path, do a save of the new path.
* Once the MATLABPATH is correct, you should be able to type commands like plotspec will now work. Also try spec


# Sampling & Aliasing


The objective of our work is to understand how system process the analog signals and storedalin discrete manner.
The aim of our objective is to discretize the signal in such a manner so that the reconstruction signal keeps all the information so that we can recover it at the end.

To check on this, Two phenomenon are very important sampling and aliasing.
One cause another, So we can concentrate on Sampling at first. We will discuss Two types of Sampling : (1) Temporal(Time Dependent)Ex: Speech signal-1D (2) Spatial (Space Dependent) Ex: Image-2D

**Theory**

Speech Signal: For experience the sampling and aliasing effect, it is better to start our investigation on Chirp signal

A linear-FM (chirp) signal is an ideal test signal to explore the concept of aliasing due to sampling. By visualizing the spectrogram of a synthesized chirp and listening to the sound, we experience the fact that a D-to-C converter cannot create output signals with frequencies higher than one half of the sampling frequency

In the demo, the chirp rate is α=1000 Hz/s, and the time duration is T=3s, so the instantaneous frequency (fi(t)=αt) goes from 1000Hz at t=0 to 3α=4000Hz at t=4s. Since the highest frequency in x(t) is 10KHz, there should be no aliasing if fs≥10KHz.
** change the FS to 10000***
The demo shows what happens for three different sampling frequencies: 10000, 5000, 2000 Hz. The first  case exhibit no aliasing and the sound rises from 1000 to 4000 Hz. In the last two cases, the sound goes up and down because fs/2<4000Hz. When fs=10000, the maximum output frequency is 5000Hz, for fs=5000, the maximum output frequency is 2500Hz, and for fs=2000, the maximum output frequency is 1000Hz.



A periodic signal is known to have a Fourier Series,which is usually described as a harmonic line spectrum because the only frequencies present in the spectrum are integer multiples of the fundamental frequency. 
In our experiment with the spectrogram, we exhibit this harmonic line characteristic. 
First we generate a periodic signal of period 10msec using a sampling rate 10000Hz. Keeping the duration 3sec and observing the 5 periods of the generated signal we can see that it is a triangular periodic signal.
We then observe the spectrum of this signal and was able to spot the harmonics. However the higher harmonics are not clearly visible in linear scale of spectrogram as they have very small amplitudes.Hence we use a logarithmic or dB scale to spot all the harmonics in spectrogram.
By changing the period from 20msec to 4msec we notice that when period is shorter the frequency separation of the harmonic lines is greater. 
# How to run Exercise
* Save the main.m
* Verify "plotspec" is accessible using toolbox.
* Run the code or publish for better visualization.
