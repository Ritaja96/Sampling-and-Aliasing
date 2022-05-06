%% Project_2
% ECE 6530 Matlab Project2 -- Ekata Mitra
% Submitted by Ekata Mitra (U1303742) in colaboration with Joshua Jacob, Elias Skinner,Ritaja Das

%% (a) First Part: a signal that “chirps” up to a very high frequency, 
% and the instantaneous frequency goes past half the sampling rate:
% A linear-FM (chirp) signal is an ideal test signal to explore 
% the concept of aliasing due to sampling. 
% By visualizing the spectrogram of a synthesized chirp 
% we experience the fact that a D-to-C converter cannot create 
% output signals:
% with frequencies higher than one half of the sampling frequency.
% Clearly, it can be said that:
% f_s/2 <f_max : case of aliasing
% If a signal goes past the half the sampling rate, 
% the frequencies that are recorded in the discrete
% approximation of the continuous signal are those of lower frequencies. 
% Since f_s/2 is the Nyquist frequency and
% the analog signal has a frequency higher than the Nyquist frequency,
% aliasing takes place.
%%  2.1
% The starting frequency is 1000 Hz at t=0s. 
% f_st = 1000 Hz. 
% The final frequency is 11,000 Hz at t = 4s.
% chirp rate = (11000-1000)/4 = 2500Hz/s.
% The instantaneous frequency (fi(t)) goes from 1000Hz at t=0s to 11000Hz
% at t=4s
% slope of instantaneous frequency = 2*mu
% formula:
% 2*mu*4 + 1000 = 11000  
% mu = (11000-1000)/8
% Given sampling frequency(f_s) = 4000Hz,
% Highest Frequency of Input Signal(f_max) : 11000Hz
%% 2.1 (a) Second Part: Determine the parameters for the signal
% Sampling Frequecy
f_s = 4000; 
% Starting Frequency
f_st = 1000; 
% Ending Frequency
f_end = 11000;
% Chirp signal Start Time
t_st = 0; 
% Chirp signal End Time
t_end = 4; 
% formula for finding Slope of instantaneous frequency:
% 2*mu = (f_end - f_st)/(t_end -t_st)
% slope =2*mu
mu = (f_end-f_st)/(2*(t_end -t_st)); 
% Time Step
dt = 1/f_s;
% Time Vector
tt = t_st:dt:t_end; 
% Random Initial Phase
phi = 2*pi*rand; 
% Code for generating chirp signal
% Signal Phase for all times
psi = 2 * pi * mu * tt.^2 + 2 * pi * f_st * tt + phi; 
% Chirp Signal
cc = real(7.7 * exp(1j*psi)); 
% plot for visualization
figure;
% function plotspec
% him = plotspec(xx,fsamp,Lsect)
%       him = handle to the image object
%        xx = input signal
%     fsamp = sampling rate
%     Lsect = section length (integer number of samples, should be power of 2)
%               amount of data to Fourier analyze at one time
% L_sect =512;
plotspec(cc, f_s, 512);
xlabel('Time(sec)');
ylabel('Frequency(Hz)');
%% 2.1(b) Generate the signal and plot in spectrogram with short section length
% The section duration is equal to (Section length / Sampling Frequency).
% Here is the spectrogram that is plotted for Lsect = 120.
% Let, L_sect = 128
% T_sect = L_sect/f_s = 128/4000s = 0.0320 s
% function plotspec
% him = plotspec(xx,fsamp,Lsect)
%       him = handle to the image object
%        xx = input signal
%     fsamp = sampling rate
%     Lsect = section length (integer number of samples, should be power of 2)
%               amount of data to Fourier analyze at one time
% L_sect =128
%%
% Code for generating chirp signal
L_sect = 128;
% Time Step
dt = 1/f_s; 
% Time Vector
tt = t_st:dt:t_end; 
% Random Initial Phase
phi = 2*pi*rand; 
% Signal Phase for all times
psi = 2 * pi * mu * tt.^2 + 2 * pi * f_st * tt + phi; 
% Chirp Signal
cc = real(7.7 * exp(1j*psi)); 

% plot for visualization
figure;

plotspec(cc, f_s, L_sect);
xlabel('Time(sec)');
ylabel('Frequency(Hz)');
%% 2.1 c
% This is because of Aliasing. 
% Here f_s/2 = 2000 Hz. 
% Therefore, 2000 Hz is the maximum frequency that can be recorded in this
% discrete approximation of the signal. At higher values,i.e.,
% when the signal frequency goes over 2000 Hz, the frequency wraps around 
% the 2000 Hz limit and decreases to 0 when f_s = 4000 Hz.
% Again, it now continues to increase to 2000 Hz,
% until the instantaneous frequency is 6000 Hz, and 
% then goes back to 0 at 8000 Hz and so on. This phenomenon is called aliasing
%% Plotspec Function for visualization
function him = plotspec(xx,fsamp,Lsect)
%PLOTSPEC   plot a Spectrogram as an image
%         (take care of all the default settings)
%  usage:   him = plotspec(xx,fsamp,Lsect)
%      him = handle to the image object
%       xx = input signal
%    fsamp = sampling rate
%    Lsect = section length (integer number of samples, should be power of 2)
%              amount of data to Fourier analyze at one time
%

% if( exist('specgram'))
% 	disp(' ')
% 	disp('??? Why are you using this function, you seem to have SPECGRAM');
% 	disp(' ')
% end
if( nargin<3 )
	Lsect = 256;
end
if( nargin<2 )
	disp('PLOTSPEC: Sampling Frequency defaulting to 8000 Hz')
	fsamp = 8000;
end
if( length(xx)<1000 )
	warning('PLOTSPEC: Signal length must be greater than 1000 to get a reasonable spectrogram')
end
Lfft = Lsect;
Noverlap = round(Lsect/2);  %-- overlap defaults to 50%
[B,F,T] = spectgr(xx,Lfft,fsamp,Lsect,Noverlap);
him = imagesc(T,F,abs(B));
axis xy
colormap(1-gray)   %-- use colormap(jet) if you like bright colors !
end