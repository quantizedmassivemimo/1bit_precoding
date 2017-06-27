 =======================================================
Simulator for "Quantized Precoding for Massive MU-MIMO"
-------------------------------------------------------
(c) 2017 Christoph Studer and Sven Jacobsson; e-mail: studer@cornell.edu and sven.jacobsson@ericsson.com  
=======================================================

# Important information:

If you are using the simulator (or parts of it) for a publication, then you MUST cite our paper:

S. Jacobsson, G. Durisi, M. Coldrey, T. Goldstein, and C. Studer, Quantized precoding for massive MU-MIMO, IEEE Trans. Commun., Jun. 2017, to appear.

and clearly mention this in your paper. 

More information about our research can be found at: https://sites.google.com/site/durisi and at http://vip.ece.cornell.edu

# How to start a simulation:

For large MU-MIMO wireless systems:

Simply run  

>> precoder_sim

which starts a simulation in a 128 BS antenna, 16 UE, massive MU-MIMO system with QPSK modulation using various precoders. The simulator returns a plot of the BER as a function of the SNR.

The following precoders are currently supported by the simulator.
  - MRTi: Maximal-ratio transmission (infinite resolution)
  - MRT: Maximal-ratio transmission (quantized)
  - ZFi: Zero-forcing (infinite resolution)
  - ZF: Zero-forcing (quantized)
  - WFi: Wiener filter (infinite resolution)
  - WF: Wiener filter (quantized)
  - EXS: Exhaustive search
  - SP: sphere precoder
  - SDR: semidefinite relaxation (rank-one approximation and randomization)
  - SQUID: squared infinity-norm relaxation with Douglas-Rachford splitting

The precoders MRT, ZF, and WF can be used together with multi-bit DACs, whereas the precoders EXS, SP, SDR, and SQUID are supported for 1-bit DACs only. It is recommended to use EXS, SP, and SDR only for smaller systems as the complexity of these precoders scales unfavorably with the number of BS antennas.  

To use the SDR precoder, you need to install CVX, which can be found here:

http://cvxr.com/cvx/download/

The simulator runs with a set of predefined simulation parameters. You can provide your own system and simulation parameters by passing your own "par"-structure (see the simulator for an example).

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.
