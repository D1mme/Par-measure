"""
# Author:    Dimme de Groot
# Date:      July 2025
# Description: This module implements the perceptual distortion measure as described in [1]. 
#              The distortion measure is referred to as the Par-measure, named after the first author. 
#              Additional details can be found in [2, 3, 4, 5].
# Important note: I obtained this Julia version from the Matlab version that can be found at https://github.com/D1mme/Par-measure . 
#
# Usage: 
#   Initialize the class as:
#       par_measure = ParMeasure(Fs, Tframe, x_ref, x_dB_ref, F_cal, Ng)
#   Compute the masking curve and weighting curves as:       
#       maskcurve, maskcurve_spl, p_par = comp_maskcurve(par_measure, masker)
#   Optionally, set all values of p_par below a certain frequency artificially high:
#       maskcurve, maskcurve_spl, p_par = comp_maskcurve(par_measure, masker; flag_low_freq=false, threshold)
#   Plot the masking curve:
#       plot_maskcurve(par_measure, masker)
#   Plot the masking curve and the disturbance:
#       plot_maskcurve(par_measure, masker, disturbance)
#   Convert between physical and digital amplitude representations:
#       amplitude_digital = physical_to_digital(par_measure, amplitude_physical)
#       amplitude_physical = digital_to_physical(par_measure, amplitude_digital)
#       Note that these amplitudes assume a straightforward mapping from the instantaneous amplitude to the instantaneous amplitude in dB SPL. 
#       In any practical scenario, this is probably not the case. But it is likely close enough
#       
# Inputs:
#   Fs:             [Hz], sample frequency              
#   Tframe:         [s], duration of the zero-padded segmentation window (typically 20 to 40 ms)
#   x_ref:          [-], reference digital value for known sound pressure level (e.g., x_ref = 1)
#   x_dB_ref:       [dB SPL], corresponding sound pressure level (e.g., x_dB_ref = 70)
#   F_cal:          [Hz], calibration frequency (typically 1000 Hz)
#   Ng:             [-], number of gammatone filters used (typically 32 to 64)
#   masker:         Discrete time-domain signal (masker) of length Nframe 
#   disturbance:    Discrete time-domain signal (disturbance) of length Nframe, optional
#   threshold:      Threshold [Hz] below which p_par values are set artificially high
#
# Outputs:
#   maskcurve:      [-], double-sided masking curve (useful for internal computations)
#   maskcurve_spl:  [dB SPL], single-sided masking curve in dB SPL
#   p_par:          [-], 1/maskcurve. To compute the Par measure, use ||diag(p_par)*epsilon||_2^2, 
#                        where epsilon is the frequency-domain disturbance  
#
# References:
# [1] Van de Par et al. A perceptual model for sinusoidal audio coding based on spectral integration, 2005. https://doi.org/10.1155/ASP.2005.1292
# [2] de Groot. A heuristic approach to spatial audio using consumer loudspeaker systems (appendix D, E, F), 2023. http://resolver.tudelft.nl/uuid:ee669571-b15a-4f0d-99cd-4c6c35acbd45
# [3] C. Taal. Prediction and Optimization of Speech Intelligibility in Adverse Conditions. PhD thesis. https://doi.org/10.4233/uuid:407c0a83-25c1-4cb3-ab4c-5a691d19f949
# [4] B. R. Glasberg and B. C. Moore, "Derivation of auditory filter shapes from notched-noise data,” Hearing Research, vol. 47, no. 1, pp. 103–138, 1990.
# [5] N. de Koeijer. Sound Zones with a Cost Function based on Human Hearing (appendix A), 2021. http://resolver.tudelft.nl/uuid:4e66fb60-2350-4010-867a-53c08a0177a7
# [6] G. Charestan et al. A Gammatone-based Psychoacoustical Modeling Approach for Speech and Audio Coding, 2001. https://www.researchgate.net/publication/2884296_A_Gammatone-based_Psychoacoustical_Modeling_Approach_for_Speech_and_Audio_Coding
#
#
# Remarks:
#   Some parts of the code could be improved, such as the factor of two in the masking curve computation. 
#   This might be resolved by refining the calibration procedure and tracking FFT normalizations more carefully. 
#   Despite these issues, the implementation aligns well with the plots in [1] and provides useful predictions for distortion audibility.
"""

using FFTW
using Plots
using SpecialFunctions   # for factorial
using Printf             # for formatted printing
using LinearAlgebra
using Random
using Statistics


mutable struct ParMeasure
    Fs::Float64
    Tframe::Float64
    x_ref::Float64
    x_dB_ref::Float64
    F_cal::Float64
    Ng::Int

    Nframe::Int
    freq_ax::Vector{Float64}
    alpha_ov_p0::Float64
    h_hat_om::Vector{Float64}
    h_hat_gamma::Matrix{Float64}
    c::Vector{Float64}

    function ParMeasure(Fs=48000.0, Tframe=0.4, x_ref=1.0, x_dB_ref=70.0, F_cal=1000.0, Ng=64)
        # Calculate nbr of samples/frame
        Nframe = round(Int, Fs*Tframe/2)*2
        println("The distortion measure expects frames of $Nframe samples")

        # Compute freq. ax in [0, Fs/2]
        freq_ax = (0:Nframe÷2) .* Fs / Nframe

        # Needed for conversion from digital to SPL [2, Appendix F], [5, Appendix A]
        alpha_ov_p0 = 10^((x_dB_ref - 20*log10(x_ref))/20)
        
        # The outer middle ear filter of Par. [3, p. 23]. Note that this is just the threshold in quiet
        h_hat_om = _compute_h_hat_om(freq_ax, alpha_ov_p0)
        
        # The Gammatone filters [1], [4]
        h_hat_gamma = _compute_h_hat_gamma(freq_ax, Ng, Fs)

        # Calibration constants [1], [6]. 
        # Note: if i recall correctly, the proof in [6] is incorrect: the calibration might not converge if you pick the wrong frequency
        # Namely, the lim C1-->infty F(C1) after eq. 6 misses a -1: thus, F(C1-->infty) could be smaller than zero, in which case the biseciton method does not converge 
        c = _compute_c(Nframe, freq_ax, h_hat_om, h_hat_gamma, F_cal, alpha_ov_p0)

        new(Fs, Tframe, x_ref, x_dB_ref, F_cal, Ng, Nframe, freq_ax,
            alpha_ov_p0, h_hat_om, h_hat_gamma, c)
    end
end


function _compute_h_hat_om(freq_ax, alpha_ov_p0)
    f = abs.(freq_ax)
    Tq = 3.64 .* (f ./ 1000) .^ (-0.8) .- 6.5 .* exp.(-0.6 .* (f ./ 1000 .- 3.3).^2) .+ 1e-3 .* (f ./ 1000).^4
    Tq .-= 20 * log10(alpha_ov_p0)
    return 10 .^ (-Tq ./ 20)
end


function _compute_h_hat_gamma(freq_ax, Ng, Fs)
    nu = 4
    kappa = 2^(nu - 1) * factorial(nu - 1) / (π * _doublefactorial(2 * nu - 3))
    fc = _erb_spacing(Fs, Ng)
    f = abs.(freq_ax)

    h_hat_gamma = zeros(Ng, length(freq_ax))
    ERB(f) = 24.7 * (4.37 * f / 1000 + 1)

    for ng in 1:Ng
        h_hat_gamma[ng, :] .= (1 .+ ((f .- fc[ng]) ./ (kappa * ERB(fc[ng]))).^2) .^ (-nu/2)
    end

    return h_hat_gamma
end


function _erb_spacing(Fs, Ng)
    E_0 = 0.0
    E_Fs = 21.4 * log10(4.37 * Fs / 2000 + 1)
    E = range(E_0, E_Fs, length=Ng)
    return 1000 / 4.37 .* (10 .^(E ./ 21.4) .- 1)
end


function _doublefactorial(n)
    p = 1
    for i in (iseven(n) ? 2 : 1):2:n
        p *= i
    end
    return p
end


function _compute_c(Nframe, freq_ax, h_hat_om, h_hat_gamma, F_cal, alpha_ov_p0)
    A_52 = 10^(52/20)/alpha_ov_p0
    A_70 = 10^(70/20)/alpha_ov_p0

    _, idx = findmin(abs.(freq_ax .- F_cal))
    println("Calibration frequency offset: $(abs(freq_ax[idx] - F_cal)) Hz")

    h_gamma = h_hat_gamma[:, idx]
    h_om = h_hat_om[idx]
    Term1 = sum(h_gamma.^2)
    Term2 = A_70^2 * h_om^2 / Nframe
    Term3 = A_52^2 * h_om^2

    c2_l = 0.0
    c2_r = 10.0
    tol = 1e-15
    f = c -> c * sum(Term3 * h_gamma.^2 ./ (Term2 * h_gamma.^2 .+ c * Term1)) - 1

    if f(c2_l)*f(c2_r) > 0
        error("Calibration failed: try increasing c2_r or changing F_cal.")
    end

    while abs(c2_r - c2_l) > tol
        c2_m = (c2_l + c2_r)/2
        if f(c2_m)*f(c2_l) < 0
            c2_r = c2_m
        else
            c2_l = c2_m
        end
    end

    c2 = (c2_l + c2_r)/2
    c1 = c2 / Nframe * Term1
    return [c1, c2]
end


function comp_maskcurve(p::ParMeasure, x::Vector{Float64}; flag_low_freq::Bool=true, threshold::Float64=30.0)
    x = collect(x)
    x_hat = fft(x) / length(x)
    x_hat = 2 .* x_hat[1:(p.Nframe ÷ 2 + 1)]  # single-sided

    tmp = zeros(length(x_hat))
    for i in 1:p.Ng
        num = p.h_hat_om.^2 .* p.h_hat_gamma[i, :].^2
        den = norm(x_hat .* p.h_hat_om .* p.h_hat_gamma[i, :])^2 / p.Nframe + p.Nframe * p.c[1]
        tmp .+= num ./ den
    end

    gsqr_SS = p.c[2] .* tmp
    gsqr_DS = 2 .* vcat(gsqr_SS, reverse(gsqr_SS[2:end-1]))

    p_par = sqrt.(gsqr_DS) ./ p.Nframe
    maskcurve = 1.0 ./ sqrt.(gsqr_DS)
    maskcurve_spl = 20 .* log10.(1.0 ./ sqrt.(gsqr_SS)) .+ 20 .* log10(p.alpha_ov_p0)

    if !flag_low_freq
        ind_thresh = findfirst(x -> x > threshold, p.freq_ax)
        if ind_thresh !== nothing
            maxval = maximum(p_par)
            p_par[1:ind_thresh] .= 100 * maxval
            p_par[end-ind_thresh+2:end] .= 100 * maxval
        end
    end

    return maskcurve, maskcurve_spl, p_par
end


function physical_to_digital(p::ParMeasure, A_SPL::Float64)
    return 10^(A_SPL / 20) / p.alpha_ov_p0
end


function digital_to_physical(p::ParMeasure, A_digital::Float64)
    return 20 * (log10(abs(A_digital)) + log10(p.alpha_ov_p0))
end


function plot_maskcurve(p::ParMeasure, masker::Vector{Float64}, disturbance::Union{Nothing, Vector{Float64}}=nothing)
    f = p.freq_ax
    Tq = 3.64 .* (f ./ 1000).^(-0.8) .- 6.5 .* exp.(-0.6 .* (f ./ 1000 .- 3.3).^2) .+ 1e-3 .* (f ./ 1000).^4

    _, maskcurve_spl, _ = comp_maskcurve(p, masker)
    masker_fft = 2 .* abs.(fft(masker)) ./ p.Nframe
    masker_fft = masker_fft[1:length(maskcurve_spl)]
    mask_SPL = 20 .* log10.(masker_fft .* p.alpha_ov_p0)

    plt = plot(f, Tq, label="Threshold in quiet", lw=2, xaxis=:log, xlabel="Frequency [Hz]", ylabel="Magnitude [dB SPL]", yticks=-10:10:100)
    plot!(f, maskcurve_spl, label="Masking curve", lw=2, linestyle=:dash)
    plot!(f, mask_SPL, label="Masker")

    if disturbance !== nothing
        disturbance_fft = 2 .* abs.(fft(disturbance)) ./ p.Nframe
        disturbance_fft = disturbance_fft[1:length(maskcurve_spl)]
        dist_SPL = 20 .* log10.(disturbance_fft .* p.alpha_ov_p0)
        plot!(f, dist_SPL, label="Disturbance")
    end

    ylims!(-10, 100)
    xlims!(20, 16000)
    
    return plt
end


function test_Par_measure()
    # Define settings
    Fs = 48000              # [Hz] Sampling frequency
    Tframe = 0.4            # [s] Frame duration
    x_ref = 1.0             # [-] Digital reference amplitude
    x_dB_ref = 70.0         # [dB SPL] Physical reference level
    F_cal = 1000.0          # [Hz] Calibration frequency
    Ng = 64                 # Number of gammatone filters

    # Create the ParMeasure object
    Par_measure = ParMeasure(Fs, Tframe, x_ref, x_dB_ref, F_cal, Ng)

    # Define example signals
    A70 = physical_to_digital(Par_measure, 70.0)
    A52 = physical_to_digital(Par_measure, 52.0)
    A50 = physical_to_digital(Par_measure, 50.0)

    Nframe = Par_measure.Nframe
    t = range(0, step=1/Fs, length=Nframe)

    masker70_1000 = A70 .* sin.(2π*1000 .* t)
    masker50_1000 = A50 .* sin.(2π*1000 .* t)
    masker50_1200 = A50 .* sin.(2π*1200 .* t)

    masker0 = zeros(length(t))
    disturbance52_1000 = A52 .* sin.(2π*1000 .* t)
    disturbance52_1200 = A52 .* sin.(2π*1200 .* t)

    masker_noise = randn(length(t))
    masker_noise = masker_noise*physical_to_digital(Par_measure, 30.0)*sqrt(Fs) # this is not entirely correct, but a quick test. 

    # Test case: 52 dB SPL disturbance and 70 dB SPL masker
    maskcurve_70, maskcurve_spl70, p_par70 = comp_maskcurve(Par_measure, masker70_1000)

    println()
    println("For a 52 dB SPL 1 kHz disturbance and a 70 dB SPL 1 kHz masker, the Par measure evaluates to: ",
        norm(p_par70 .* fft(disturbance52_1000))^2, " (should be about one)")

    println("The maskcurve should, up to a normalization, be equal to p_par. Validation: ",
        norm(p_par70 * Nframe .- 1.0 ./ maskcurve_70), " (should be about zero)")

    # Plotting masking curves
    plt = plot_maskcurve(Par_measure, masker70_1000)
    title!("Predicted masking curve for a 70 dB SPL masker (1000 Hz)", titlefont = 12)
    display(plt)

    plt = plot_maskcurve(Par_measure, masker50_1000)
    title!("Predicted masking curve for a 50 dB SPL masker (1000 Hz)", titlefont = 12)
    display(plt)

    plt = plot_maskcurve(Par_measure, masker50_1200)
    title!("Predicted masking curve for a 50 dB SPL masker (1200 Hz)", titlefont = 12)
    display(plt)

    plt = plot_maskcurve(Par_measure, masker0)
    title!("Predicted masking curve when no masker is present.", titlefont = 12)
    display(plt)

    plt = plot_maskcurve(Par_measure, masker70_1000, disturbance52_1000)
    title!("Predicted masking curve with 70 dB SPL masker (1000 Hz)\nJust (un)noticeable disturbance.", titlefont = 12)
    display(plt)

    plt = plot_maskcurve(Par_measure, masker70_1000, disturbance52_1200)
    title!("Predicted masking curve with 70 dB SPL masker (1000 Hz)\nAudible disturbance.", titlefont = 12)
    display(plt)

    plt = plot_maskcurve(Par_measure, masker_noise)
    title!("Predicted masking curve with noise as masker.", titlefont = 12)
    display(plt)
end

