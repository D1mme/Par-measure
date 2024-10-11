%Author:    Dimme de Groot
%Date:      April 2024
%Descr:     Simple example of how to use Par-measure (see [1] inside optimisation problem (CVX, see [2])
%           This function attempt to increase the loudness of signals while decreasing their maximum amplitude.
%           This is based on the work presented in [3].
%Sources:   
%   [1] Van de Par et al. A perceptual model for sinusoidal audio coding based on spectral integration, 2005. https://doi.org/10.1155/ASP.2005.1292
%   [2] CVX Research, Inc. CVX: Matlab software for disciplined convex programming, version 2.2, Build 1184. URL: http://cvxr.com/cvx
%   [3] Jeannerot et al. Increasing Loudness in Audio Signals: A perceptually motivated approach to preserve audio quality, IEEE ICASSP 2022. 

clear all
close all

%user setting: maximum distortion dPar, original audio file and window length 
dPar = 40;                       %Maximum allowable distortion
example = "Example_audio_1";
audiofile = "Data/" + example + "/reference.wav"; %Reference audio file
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise problem setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read reference audio
[s_ref, Fs] = audioread(audiofile); 

%Initialise par measure
Twin = length(s_ref)/Fs;     
x_ref = 1; x_dB_ref = 90;   %[-], [dB SPL]; the reference value in digital and physical domain
F_cal = 400;               %[Hz], the calibration frequency. Note that in the report this corresponds to f_m
Ng = 64;                    %[-], the number of gammatone filters used
Par_meas = par_measure(Fs, Twin, x_ref, x_dB_ref, F_cal, Ng);   

%Compute Par_meas.Nframe x Par_meas.Nframe DFT matrix
W = dftmtx(Par_meas.Nframe);
s_ref = [s_ref; zeros(Par_meas.Nframe-length(s_ref),1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over optimisaiton problem on a frame by frame basis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for output
s_out = fnc_Jeannerot(s_ref, Par_meas, W, dPar);   %solve optimisation problem
      
%%%%%%%%%%%%%%%%
% Store result %
%%%%%%%%%%%%%%%%
audiowrite("Data/" + example + "/loudness_percep_"+num2str(dPar)+".wav", s_out/max(abs(s_out)), Fs);        


function s = fnc_Jeannerot(s_ref, par_meas, W, dPar)
    %this function computes one instance of the optimisation problem proposed by [3] for increased loudness
    [~, ~, p_par] = par_meas.comp_maskcurve(s_ref, false, 30);   %Compute masking curve
    Lwin = length(p_par);
    P_par = diag(p_par);  %Put masking curve on a diagonal matrix  

    cvx_solver SDPT3      %Set the solver to SDPT3: this is the default solver and comes with CVX!
    cvx_begin 
        variable s(Lwin)
        minimise norm(s, inf) %Minimise the infinity norm (maximum absolute value)
        subject to                                          
            norm(P_par*W*(s-s_ref)) <= sqrt(dPar);    %Subject to keeping the distortion as measured by the Par-measure limited    
    cvx_end 
end