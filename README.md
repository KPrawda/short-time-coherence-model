# Short-time Coherence Model
 
Complementary code for the publication "Short-time Coherence Between Repeated Room Impulse Response Measurements" by K. Prawda, S. J. Schlecht, and V. Välimäki, published by the Journal of the Acoustical Society of America on 13.08.2024.<br>
See the article here: https://doi.org/10.1121/10.0028172



**Abstract:**<br>
Room impulse responses (RIRs) vary over time due to fluctuations in atmospheric temperature, humidity, and pressure. This can introduce uncertainties in room transfer-function measurements, which are challenging to account for. Previous methods of identification and compensation of time variance focus on systematic atmospheric changes and do not apply to subtle discrepancies in RIRs. In this work, we address this problem by proposing a model of short-time coherence between repeated RIR measurements as an indicator of time-frequency similarity and as a measure of time-variance-induced changes in RIRs. Atmospheric changes cause fluctuation in sound speed, which, in turn, results in variation in the time-of-arrival of sound reflections following a Generalized Wiener process. We show that the short-time coherence decreases exponentially with the reflection-path length and propose volatility as a single model parameter determining the coherence decay rate. The proposed model is validated on simulations and measurements, showing applicability in indoor scenarios. The method reliably estimates volatility of $10^{-6}\textrm{s}/\sqrt{\textrm{s}}$ as measured under laboratory conditions. We exemplify the utility of short-time coherence loss by predicting the high-frequency energy loss stemming from RIR averaging. The proposed method is useful in assessing the uncertainty of RIR measurements, especially when repeated measurements are compared or averaged.

<br><br>
![Two ISM examples with different speed of sound voxels](https://github.com/KPrawda/short-time-coherence-model/blob/main/Figures/Fig1.jpg)
A two-dimensional (2-D) view of a simulated room (black contour) with a sound source S (black dot), its first-order images IS1 and IS2 (gray dots), and a sound receiver R (white dot) with the speed of sound distribution (a) c(x) and (b) c'(x). The spatial distribution of the speed of sound is marked with colored voxels. The reflection paths are marked with dashed lines.

<br><br>
![RIR differences stemming from time variance](https://github.com/KPrawda/short-time-coherence-model/blob/main/Figures/Fig2.png)
RIRs of simulated rooms from the figure above, obtained with ISM using two different c(x) and c'(x). The difference between c(x) and c'(x) leads to discrepancies in reflections that grow over RIR duration.


The codes reproduce all the figures from the paper


**Cite as:**<br>
@article{prawda2024_coherence,<br>
author = {Prawda, Karolina  and Schlecht, Sebastian J.  and Välimäki, Vesa},<br>
title = {Short-time Coherence Between Repeated Room Impulse Response Measurements},<br>
journal = {J. Acoust. Soc. Am.},<br>
volume = {156},<br>
number = {2},<br>
pages = {1017--1028},<br>
year = {2024},<br>
doi = {https://doi.org/10.1121/10.0028172}}<br>
