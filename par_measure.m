%Author:    Dimme de Groot
%Date:      May 2024
%Descr:     This class implements the perceptual distortion measure as described in [1]. I refer to this distortion measure as the Par-measure, after the first author. 
%           Some more details can be found in [2, 3, 4, 5]
%
%Usage:     Initialize the class as:
%               Par_measure = par_measure(Fs, Tframe, x_ref, x_dB_ref, F_cal, Ng)
%           Compute the masking curve and weighting curves as:       
%               [maskcurve, maskcurve_spl, p_par] = par_measure.comp_maskcurve(masker)
%           In some scenarios, it might be convenient to set all values of p_par below a certain frequency high. In that case you can use
%               [maskcurve, maskcurve_spl, p_par] = par_measure.comp_maskcurve(masker, false, threshold)
%               with threshold the frequency value below which p_par is set high
%           Plot the masking curve as:           
%               Par_measure.plot_maskcurve(masker)
%           Plot the masking curve and the disturbance as               
%               Par_measure.plot_maskcurve(masker, disturbance)
%           Compute the digital representation amplitude given a physical representation amplitude (in dB SPL) and vice versa using: 
%               amplitude_digital = Par_measure.physical_to_digital(amplitude_physical)
%               amplitude_physical = Par_measure.digital_to_physical(amplitude_digital) 
%
%Inputs:        Fs:         [Hz], sample frequency              
%               Tframe:     [s], the time of the zero-padded segmentation window                        (typically 20 to 40 ms)
%               x_ref:      [-], a reference digital value for which you know the sound pressure level  (for example, x_ref = 1)
%               x_dB_ref:   [dB SPL], the corresponding sound pressure level                            (for example, x_dB_ref = 70)
%               F_cal:      [Hz], a calibration frequency.                                              (typically 1000 Hz)
%               Ng:         [-], the number of gammatone filters used.                                  (typically 32 to 64)
%               
%               masker:         the discrete time-domain signal (masker) of length Nframe 
%               disturbance:    the discrete time-domain signal (distrubance) of length Nframe, optional! 
%               threshold:      the threshold [Hz] below which the values of p_par should be set high artificially. 
%                                          
%Outputs:       maskcurve:      [-], the double sided masking curve. I dont really know a scenario where this one is usefull to view
%               maskcurve_spl:  [dB SPL], the single sided masking curve in dB SPL
%               p_par:          [-], 1/maskcurve. To compute the par measure, you use ||diag(p_par)*epsilon||_2^2, with epsilon the freq. domain disturbance  
%
%Sources:
%[1] Van de Par et al. A perceptual model for sinusoidal audio coding based on spectral integration, 2005. https://doi.org/10.1155/ASP.2005.1292
%[2] de Groot.  A heuristic approach to spatial audio using consumer loudspeaker systems (appendix D, E, F), 2023. http://resolver.tudelft.nl/uuid:ee669571-b15a-4f0d-99cd-4c6c35acbd45
%[3] C. Taal. Prediction and Optimization of Speech Intelligibility in Adverse Conditions. PhD thesis. https://doi.org/10.4233/uuid:407c0a83-25c1-4cb3-ab4c-5a691d19f949
%[4] B. R. Glasberg and B. C. Moore, "Derivation of auditory filter shapes from notched-noise data,” Hearing Research, vol. 47, no. 1, pp. 103–138, 1990.
%[5] N. de Koeijer. Sound Zones with a Cost Function based on Human Hearing (appendix A), 2021. http://resolver.tudelft.nl/uuid:4e66fb60-2350-4010-867a-53c08a0177a7
%[6] G. Charestan et al. A Gammatone-based Psychoacoustical Modeling Approach for Speech and Audio Coding, 2001. https://www.researchgate.net/publication/2884296_A_Gammatone-based_Psychoacoustical_Modeling_Approach_for_Speech_and_Audio_Coding
%
%Remarks:
%   There are some parts of the code where I'm not really happy with. For example, I have a weird factor of two in the computation of the masking curve.  
%   This is likely all resolvable by changing the calibration procedure and tracking the FFT normalisations better. 
%   Anyway, through comparison with the plots in [1] and by listening to signals predicted to have an inaudible or only slightly audible distortion, 
%   I feel like the implementation is at least mostly correct and usefull as a predictor of distortions.

classdef par_measure
    properties
        Fs                      %[Hz], sample frequency
        Tframe                  %[s], the time of the zero-padded segmentation window
        x_ref                   %[-], the reference value in digital domain
        x_dB_ref                %[dB SPL] the reference value in physical domain
        F_cal                   %[Hz], the calibration frequency
        Ng                      %[-], the number of gammatone filters used

        Nframe                  %[-],   the number of samples per window --> obtained from Tframe, even number
        freq_ax                 %[Hz],  the frequency axis 
        alpha_ov_p0             %[-],   a constant useful in calculating the SPL in dB
        h_hat_om                %[-],   the outer middle ear filter (misleading! it is the inverse of the threshold in quiet) used in Par
        h_hat_gamma             %[-],   the gammatone filter used in Par
        c                       %[c1, c2] [-], the calibration constants used in Par
    end
        
    methods
        function obj = par_measure(Fs, Tframe, x_ref, x_dB_ref, F_cal, Ng)
            if nargin == 0
                obj.Fs = 48000;
                obj.Tframe = 0.4;              
                obj.x_ref = 1;
                obj.x_dB_ref = 70;   
                obj.F_cal = 1000;               
                obj.Ng = 64;                    
            else
                obj.Fs = Fs;
                obj.Tframe = Tframe;              
                obj.x_ref = x_ref;
                obj.x_dB_ref = x_dB_ref;   
                obj.F_cal = F_cal;               
                obj.Ng = Ng;    
            end
       
            obj.Nframe = methodNframe(obj);
            obj.freq_ax = methodFreq_ax(obj);
            obj.alpha_ov_p0 = methodAlpha_ov_p0(obj);
            obj.h_hat_om = methodH_hat_om(obj);
            obj.h_hat_gamma = methodH_hat_gamma(obj);
            obj.c = methodC(obj);
        end

        %Calculate number of samples/frame
        function Nframe = methodNframe(obj)
            Nframe = round(obj.Fs*obj.Tframe/2)*2;
            disp("The distortion measure expects frames of " + num2str(Nframe) + " samples")
        end

        %calculate frequency axis in [0, Fs/2]
        function freq_ax = methodFreq_ax(obj)
            k = 0:obj.Nframe/2;                 
            freq_ax = k*obj.Fs/obj.Nframe;      
        end

        %Needed for conversion from digital to SPL [2, Appendix F], [5, Appendix A]
        function alpha_ov_p0 = methodAlpha_ov_p0(obj)
            alpha_ov_p0 = 10^((obj.x_dB_ref - 20*log10(obj.x_ref))/20);
        end

        %The outer middle ear filter of Par. [3, p. 23]. Note that this is just the threshold in quiet
        function h_hat_om = methodH_hat_om(obj)                                     
            f = abs(obj.freq_ax);                                                           %Note: equation only holds for nonnegative 
            Tq = 3.64*(f/1000).^(-0.8)-6.5*exp(-0.6*(f/1000-3.3).^2)+10^(-3)*(f/1000).^4;   %Threshold in quit [dB SPL]
            Tq = Tq -  20*log10(obj.alpha_ov_p0);                                           %Go from dB SPL to dB
            h_hat_om = 10.^-(Tq/20);                                                        %From dB to ampltiude       
        end

        %The Gammatone filters [1], [4]
        function h_hat_gamma = methodH_hat_gamma(obj)
            nu = 4;                                             %the filter order
            kappa = 2^(nu-1)*factorial(nu-1)/(pi*dfac(2*nu-3)); %normalisation term      
            fc = spacing(obj.Fs, obj.Ng);                       %spacing between the different center frequencies of the gammatone filters
    
            f = abs(obj.freq_ax);   
            
            h_hat_gamma = zeros(obj.Ng, obj.Nframe/2+1);        %single-sided

            ERB = @(f) 24.7*(4.37*f/1000+1);                    %see [4]
            for ng=1:obj.Ng
                h_hat_gamma(ng,:) = (1+((f-fc(ng))/(kappa*ERB(fc(ng)))).^2 ).^(-nu/2);  %see [4]
            end
            
            function out = dfac(in)             %double factorial
                out = 1;
                if mod(in,2) == 0
                    for i=2:2:in
                        out = out*i;
                    end
                else
                    for i=1:2:in
                        out = out*i;
                    end
                end
            end

            function fc = spacing(Fs, Ng)
                E_0 = 0.0;                        %
                E_Fs = 21.4*log10(4.37*Fs/2000 + 1); %Fs/2
                E = linspace(E_0, E_Fs, Ng);
                fc = 1000/4.37 * (10.^(E/21.4) - 1);
            end
        end

        %Calibration constants [1], [6]. 
        %Note: if i recall correctly, the proof in [6] is incorrect: the calibration might not converge if you pick the wrong frequency
        %Namely, the lim C1-->infty F(C1) after eq. 6 misses a -1: thus, F(C1-->infty) could be smaller than zero, in which case the biseciton method does not converge 
        function c = methodC(obj)
            %find the amplitudes required for calibration
            A_52 = 10^(52/20)/obj.alpha_ov_p0;
            A_70 = 10^(70/20)/obj.alpha_ov_p0;
                
            %find the filters associated with f_1 and f_2   (which are the same in this implementation)
            [M, Indx] = min(abs(obj.freq_ax - obj.F_cal));
            disp("Calibration frequency offset: " + num2str(M) + " Hz") 
            h_hat_gamma_f1 = obj.h_hat_gamma(:,Indx);
            h_hat_om_f1 = obj.h_hat_om(Indx);

            h_hat_gamma_f2 = h_hat_gamma_f1;
            h_hat_om_f2 = h_hat_om_f1;
       
            %Some terms which can be precalculated
            Term1 = sum(h_hat_gamma_f1.^2);
            Term2 = A_70^2*h_hat_om_f2^2/obj.Nframe;
            Term3 = A_52^2*h_hat_om_f2^2; 

            %Bisection Method:
            tol = 1e-15;
            f_c2t = 2*tol;
            
            c2_l = 0;       %inital guess
            c2_r = 10;      %initial guess 

            f_c2l = fnc_f_c2(c2_l, Term1, Term2, Term3, h_hat_gamma_f2);
            f_c2r = fnc_f_c2(c2_r, Term1, Term2, Term3, h_hat_gamma_f2);
  
            if f_c2l*f_c2r > 0
                disp('ABORT: solution not in guessed interval. You can fix this by increasing c2_r, though it is more likely that a weird calibration frequency is used')
            else   
                while abs(f_c2t)>tol
                    c2_t = (c2_l+ c2_r)/2;
                    f_c2t = fnc_f_c2(c2_t, Term1, Term2, Term3, h_hat_gamma_f2) ;
                    if f_c2t*f_c2l < 0
                        c2_r = c2_t;
                        f_c2r = f_c2t;
                    else
                        c2_l = c2_t;
                        f_c2l = f_c2t;
                    end
                end
            
            end

            c2 = c2_t;
            c1 = c2/obj.Nframe * Term1;
            c = [c1, c2];
            function f_c2 = fnc_f_c2(c2, Term1, Term2, Term3, h_hat_gamma_f2) 
               f_c2 = c2*sum( Term3*h_hat_gamma_f2.^2 ./ (Term2*h_hat_gamma_f2.^2 + c2*Term1) )-1;
            end
        end

        %Functions for users of the measure
        function [maskcurve, maskcurve_spl, p_par] = comp_maskcurve(obj, x, flag_low_freq, threshold)
            if nargin == 2
                flag_low_freq = true;
            end
            if nargin == 3 && ~flag_low_freq
                threshold = 30;%Hz
                [~, INDXthresh] = min(abs(obj.freq_ax-threshold));
            end
           
            if nargin == 4 && ~flag_low_freq
                [~, INDXthresh] = min(abs(obj.freq_ax-threshold));
            end
            x = x(:).'; %ensure row
            x_hat = fft(x)/length(x);
        
            %find g^2
            x_hat = 2*x_hat(1:obj.Nframe/2+1);  %convert to single sided
            tmp = zeros(1,obj.Nframe/2+1);
            for i = 1:obj.Ng
                num = obj.h_hat_om.^2.*obj.h_hat_gamma(i,:).^2;
                den = norm(x_hat.*obj.h_hat_om.*obj.h_hat_gamma(i,:))^2/obj.Nframe+obj.Nframe*obj.c(1); 
                tmp = tmp+num./den;
            end
            gsqr_SS = obj.c(2)*tmp;                             %=g^2, single sided
            gsqr_DS = 2*[gsqr_SS, fliplr(gsqr_SS(2:end-1))];    %=g^2, double sided (times 2, huh!)
            
            p_par = sqrt(gsqr_DS)/(obj.Nframe);                 %The division by obj.Nframe is needed due to FFT stuff. It could also be part of the calibration procedure.
            maskcurve = 1./sqrt(gsqr_DS);                       %DOUBLE SIDED, [-]     
            maskcurve_spl = 20*log10(1./sqrt(gsqr_SS))+20*log10(obj.alpha_ov_p0);    %SINGLE SIDED, [dB SPL]
        
            %Artificially set weights for low frequencies high (i.e. discourage disturbances in low frequencies	
            if ~flag_low_freq
            	val = max(p_par);
                if INDXthresh ~= 0
		            p_par(1:INDXthresh) = 100*val;
                    p_par(end-INDXthresh+2:end) = 100*val;
                end
            end
        end
    
        function A_digital = physical_to_digital(obj, A_SPL) %digital <--> physical representation
            %Function to go from a physical sound level (dB SPL) to a digital representation %
            A_digital = 10^(A_SPL/20)/obj.alpha_ov_p0;
        end
        function A_physical = digital_to_physical(obj, A_digital) %digital <--> physical representation
            %Function to go from a physical sound level (dB SPL) to a digital representation %
            A_physical = 20*(log10(abs(A_digital))+log10(obj.alpha_ov_p0));
        end

        function plot_maskcurve(obj, masker, disturbance)
            if nargin == 2
                dist_flag = 0;
            else
                dist_flag = 1;
            end

            f = obj.freq_ax;        
            Tq = 3.64*(f/1000).^(-0.8)-6.5*exp(-0.6*(f/1000-3.3).^2)+10^(-3)*(f/1000).^4;   %Threshold in quit [dB SPL]
            [~, maskcurve_spl, ~] = comp_maskcurve(obj, masker);                            %masking curve [dB SPL]
            masker = 2*abs(fft(masker))/obj.Nframe;                                         %masker, single sided
            masker = masker(1:length(maskcurve_spl));                                       %"                  "
            mask_SPL = 20*log10(masker*obj.alpha_ov_p0);                                    %to dB SPL
            
            if dist_flag == 1
                disturbance = 2*abs(fft(disturbance))/obj.Nframe;                               %disturbance, single sided
                disturbance = disturbance(1:length(maskcurve_spl));                             %"                  "
                dist_SPL = 20*log10(disturbance*obj.alpha_ov_p0);                               %"       "
            
                figure
                semilogx(obj.freq_ax, Tq, 'linewidth', 2)
                hold on
                semilogx(obj.freq_ax, maskcurve_spl, 'linewidth', 2, 'linestyle', '--')
                semilogx(obj.freq_ax, mask_SPL)
                semilogx(obj.freq_ax, dist_SPL)
                grid on
                legend('Threshold in quiet', 'Masking curve', 'Masker', 'Disturbance')
                xlabel('Frequency [Hz]')
                ylabel('Magnitude [dB SPL]')
                ylim([-10 100])
                xlim([50 16000])
            else
                figure
                semilogx(obj.freq_ax, Tq, 'linewidth', 2)
                hold on
                semilogx(obj.freq_ax, maskcurve_spl, 'linewidth', 2, 'linestyle', '--')
                semilogx(obj.freq_ax, mask_SPL)
                grid on
                legend('Threshold in quiet', 'Masking curve', 'Masker')
                xlabel('Frequency [Hz]')
                ylabel('Magnitude [dB SPL]')
                ylim([-10 100])
                xlim([50 16000])

            end
        end
    end
end

