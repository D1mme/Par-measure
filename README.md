# Par-measure
In this repository you can find the MATLAB code for a perceptual distortion measure based on auditory masking. Such a measure can predict if a human can notice the difference between two acoustic signals. 
In particular, this code reflects my interpretation of the distortion measure proposed by Van de Par et al. in [1]. I refer to it as the Par-measure. The major strength of the Par-measure is that, once the masker is fixed, it can be expressed as a weighted $l_2$ norm of the disturbance. This allows for incorporating the measure in optimisation problems in which it is important to keep the compuational-complexity limited.

Concretely, consider an actual acoustic signal $s$ and a reference acoustic signal $s^\star$. The difference between the two signals is $\epsilon=s-s^\star$. The Par-measure computes the distortion as

$$ D(s, \epsilon) = d. $$

If $d \leq 1$, it is predicted that a human observer can not notice the difference between $s$ and $s^\star=s-\epsilon$ in the presence of $s$. The acoustic signal $s$, which limits the audibility of the distortion, is called the masker. The Par-measure operates in the frequency domain on short-time frames of about 20 to 40 ms [1].

This repository includes:
- A report which explains the background of the Par-measure, the functionality of the code and the provided examples;
- The MATLAB code ``par_measure.m`` which can be used to compute the masking curves;
- The MATLAB code ``example1_basics.m``, which shows the functionality of the code;
- The MATLAB code ``example2a_loudness_increase.m``, which uses the Par-measure to increase the loudness of acoustic signals while keeping the perceived distortion limited. This example is directly based on the work done by Jeannerot et al. in [2] and, in my opinion, is a nice example of how to use the Par-measure inside an optimisation framework;
- The MATLAB code ``example2b_hard_clipping.m``, which uses hard clipping to perform the loudness icnrease of acoustic signals;
- The .wav files corresponding to `example_2a_loudness_increase.m` and `example_2b_hard_clipping.m` can be found in the `Data` folder. As reference files I used two different signals. The results of the first signal can be found in `Data/Example_audio_1` and the results for the second signal are found in `Data/Example_audio_2`. The files `loudness_hard_<>.wav` correspond to hard clipping with parameter `<>`. The files `loudness_percep_<>.wav` correspond to using the Par-measure with parameter $d$ equal to `<>`. 
	- I do not know an easy way to compare the two types of parameters, so the best approach is to listen to a few audio files and try to compare the perceived distortion and the perceived loudness increase.
	- Note that the maximum amplitude of the audio in each file is equal! 

## Contact
Feel free to leave a message in case you have questions, find mistakes, or have any other comments! 

## Dependencies
The code was tested on MATLAB R2023b and MATLAB R2024b on Ubuntu 23.10. I think the Par-measure can be used using default MATLAB functionality, but if I notice that certain packages are needed I will list them here. 

To run Example 2, you need [CVX](https://cvxr.com/cvx/). It should be straightforward to rewrite the examples to use other optimisation packages such as [YALMIP](https://yalmip.github.io/) or [CVXPY](https://www.cvxpy.org/), though the latter requires interfacing with Python.

## Keywords
Acoustic signal processing, perceptual distortion, distortion measure, Par-measure, psychoacoustical modelling, auditory masking

## Note
There is also a Python implementation by Niels de Koeijer (Bang & Olufsen)!! See [link](https://github.com/nielsdekoeijer/libdetectability). This implementation can be used in PyTorch as a loss function. 

## Bibiliography
[1] van de Par, S., Kohlrausch, A., Heusdens, R. et al. A Perceptual Model for Sinusoidal Audio Coding Based on Spectral Integration. EURASIP J. Adv. Signal Process. 2005, 317529 (2005). https://doi.org/10.1155/ASP.2005.1292

[2] A. Jeannerot, N. de Koeijer, P. Martínez-Nuevo, M. B. Møller, J. Dyreby and P. Prandoni, "Increasing Loudness in Audio Signals: A Perceptually Motivated Approach to Preserve Audio Quality," ICASSP 2022 - 2022 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), Singapore, Singapore, 2022, pp. 1001-1005. https://doi.org/10.1109/ICASSP43922.2022.9747589. Alternative [ARXIV](https://arxiv.org/abs/2202.08183) link.
