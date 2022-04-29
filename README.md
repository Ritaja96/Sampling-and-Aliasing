# Sampling & Aliasing


The objective of our work is to understand how system process the analog signals and storedalin discrete manner.
The aim of our objective is to discretize the signal in such a manner so that the reconstruction signal keeps all the information so that we can recover it at the end.

To check on this, Two phenomenon are very important sampling and aliasing.
One cause another, So we can concentrate on Sampling at first. We will discuss Two types of Sampling : (1) Temporal(Time Dependent)Ex: Speech signal-1D (2) Spatial (Space Dependent) Ex: Image-2D

Sampling Theory

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
