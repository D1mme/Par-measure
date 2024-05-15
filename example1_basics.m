%Author:    Dimme de Groot
%Date:      19-03-2024
%Descr:     Example of implementation and functionality of Par measure object
%Sources:
%   [1] Van de Par et al. A perceptual model for sinusoidal audio coding based on spectral integration, 2005. https://doi.org/10.1155/ASP.2005.1292

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define settings used in Par; calibrate Par %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 48000;                 %[Hz],  Sampling frequency
Tframe = 0.04;              %[s],   the time of the input frames
x_ref = 1; x_dB_ref = 65;   %[-],[dB SPL]; the reference value in digital and physical domain
F_cal = 1000;               %[Hz],  The calibration frequency. 
Ng = 64;                    %[-],   The number of gammatone filters used

%Create the object. This prints: 
%	(a) how much the actual calibration deviates from F_cal (this can safely be ignored) and (b) the frame length (in samples) of the input frames  
par_measure = Par_measure(Fs, Tframe, x_ref, x_dB_ref, F_cal, Ng);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define some example signals %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define digital amplitudes relating to 70 dB SPL, 52 dB SPL, and 50 dB SPL
A70 = par_measure.physical_to_digital(70); %Amplitude for digital signal
A52 = par_measure.physical_to_digital(52); %Idem
A50 = par_measure.physical_to_digital(50); %Idem

% Compute sinusoids with given amplitudes and a frequency of 1 kHz
t = 0:1/Fs:par_measure.Nframe/Fs - 1/Fs;    %[s], time axis
masker70_1000 = A70*sin(2*pi*1000*t);       %[-], 70 dB SPL sinusoid of 1 kHz
masker50_1000 = A50*sin(2*pi*1000*t);       %[-], 50 dB SPL sinusoid of 1 kHz
masker50_1200 = A50*sin(2*pi*1200*t);       %[-], 52 dB SPL sinusoid of 1.2 kHz

masker0 = 0*sin(2*pi*1000*t);               %[-], -infty dB SPL sinusoid (i.e. all zeros)
disturbance52_1000 = A52*sin(2*pi*1000*t);  %[-], 52 dB SPL sinusoid of 1000 Hz
disturbance52_1200 = A52*sin(2*pi*1200*t);  %[-], 52 dB SPL sinusoid of 1200 Hz


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify that for a 52 dB SPL disturbance and a 70 dB SPL masker, the Par measure evaluates to 1.                   %
% The measure should be calibrated for that. See the original paper [1]. This disturbance should be just noticeable %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[maskcurve_70, maskcurve_spl70, p_par70] = par_measure.comp_maskcurve(masker70_1000); %This function gives the double sided maskcurve, 1/maskcurve (p_par), and single-sided maskcurve in dB SPL
disp(" ")
disp("For a 52 dB SPL 1 kHz disturbance and a 70 dB SPL 1 kHz masker, the Par measure evaluates to: " + num2str(norm(p_par70.*fft(disturbance52_1000))^2) + "; The result should be about one")
disp("The maskcurve should, up to a normalisation, be equal to p_par. Validation: " + norm(p_par70*par_measure.Nframe-1./maskcurve_70) + "; The result should be about zero!")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below some plots are given for different maskers (50 dB SPL, 70 dB SPL, no masker %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot predicted masking curve (70 dB SPL)
par_measure.plot_maskcurve(masker70_1000)
title('Predicted masking curve for a 70 dB SPL masker (1000 Hz).')

%Plot predicted masking curve (50 dB SPL)
par_measure.plot_maskcurve(masker50_1000)
title('Predicted masking curve for a 50 dB SPL masker (1000 Hz).')

%Plot predicted masking curve (50 dB SPL)
par_measure.plot_maskcurve(masker50_1200)
title('Predicted masking curve for a 50 dB SPL masker (1200 Hz).')

%Plot predicted masking curve when no masker is present. This masking curve should equal the threshold in quiet.
par_measure.plot_maskcurve(masker0)
title('Predicted masking curve when no masker is present.')

%Plot predicted masking curve and disturbance. For this example, the disturbane should lie on top of the masking curve.
par_measure.plot_maskcurve(masker70_1000, disturbance52_1000)
title({'Predicted masking curve for a 70 dB SPL masker (1000 Hz).'...
    'Just (un)noticeable disturbance.'})

%Plot predicted masking curve and disturbance. For this example, the disturbance should lie above the masking curve,
par_measure.plot_maskcurve(masker70_1000, disturbance52_1200)
title({'Predicted masking curve for a 70 dB SPL masker (1000 Hz).'...
    'Audible disturbance.'})
