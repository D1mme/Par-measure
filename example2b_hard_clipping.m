%Author:    Dimme de Groot
%Date:      May 2024
%Descr:     Hard clipping of acoustic signal. 
clear all
close all

audiofile = "Data/reference.wav";
lambda = 1; %maximum ampltidue
beta = 3;   %desired increase

%Read audio file and ampltify to desired elvel
[s_ref, Fs] = audioread(audiofile); 
s_loud = s_ref*beta;

%Perform hard clipping
s_clip = min(abs(s_loud), lambda).*sign(s_loud);

%Write result
audiowrite("Data/loudness_hard_"+num2str(beta)+".wav", s_clip, Fs);