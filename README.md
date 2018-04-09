# Simulator for "Quantized Precoding for Massive MU-MIMO"
(c) 2018 Christoph Studer and Sven Jacobsson;
e-mail: studer@cornell.edu & sven.jacobsson@ericsson.com

More information about our research can be found at [http://vip.ece.cornell.edu] and [https://sites.google.com/site/durisi].

### Important information

If you are thinking of contacting us, please do not e-mail the author to ask for download instructions, installation guidelines, or the simulator itself. The code itself is well-documented and provides the essential information about the code. Note that we will NOT help to debug user-generated code that was not included in the provided software package. If, however, you notice a bug in our code, please be so kind to contact the Author.

The software package is supplied "as is," without any accompanying support services, maintenance, or future updates. We make no warranties, explicit or implicit, that the software contained in this package is free of error or that it will meet your requirements for any particular application. It should not be relied on for any purpose where incorrect results could result in loss of property, personal injury, liability or whatsoever. If you do use our software for any such purpose, it is at your own risk. The authors disclaim all liability of any kind, either direct or consequential, resulting from your use of these programs.

If you are using the simulator (or parts of it) for a publication, then you *must* cite our paper:

S. Jacobsson, G. Durisi, M. Coldrey, T. Goldstein, and C. Studer, "Quantized precoding for massive MU-MIMO", IEEE Trans. Commun., vol. 65, no. 11, pp. 4670--4684, Nov. 2017.

and clearly mention this in your paper.

### How to start a simulation:

Simply run

```sh
precoder_sim
```

which starts a simulation for a massive MU-MIMO system (with 128 BS antennas and 16 users) with QPSK modulation using MRT and ZF precoding (for 1-bit quantization and infinite resolution), as well as the SQUID algorithm proposed in the paper.

The following precoders are currently supported by the simulator:
  - BB-1: 1-bit branch-and-bound
  - MRT: maximal-ratio transmission (quantized)
  - MRT_inf: maximal-ratio transmission (infinite resolution)
  - EXS: exhaustive search
  - SDR: semidefinite relaxation (rank-one approximation and randomization)
  - SP: sphere precoding
  - SQUID: squared infinity-norm Douglas-Rachford splitting
  - WF: Wiener filter (quantized)
  - WF_inf: Wiener filter (infinite resolution)
  - ZF: zero-forcing (quantized)
  - ZF_inf: zero-forcing (infinite resolution)

To use the SDR precoder, you need to install CVX, which can be found here: http://cvxr.com/cvx/download/

The linear precoders MRT, ZF, and WF can be used together with multi-bit quantization, whereas the nonlinear precoders BB-1, EXS, SDR, SP, and SQUID are supported only for 1-bit quantization. For ill-conditioned problems, the convergence of SQUID can be improved by tuning (see detailed comments in the code). It is recommended to use BB-1, EXS, SP, and SDR only for smaller systems, as the complexity of these precoders scales unfavorably with the number of BS antennas.

If you are using any of the precoders BB-1, EXS, and SP for a publication, then you should *also* cite our paper:

S. Jacobsson, W. Xu, G. Durisi, and C. Studer, “MSE-optimal 1-bit precoding for multiuser MIMO via branch and bound,” in Proc. IEEE Int. Conf. Acoust., Speech, Signal Process. (ICASSP), Calgary, Canada, Apr. 2018, to appear.


The simulator runs with predefined parameters. You can specify your own system and simulation parameters by passing your own "par"-structure (see the simulator for an example). Note that we use default parameters for the considered system configuration; if you want to run the simulation with different parameters, then please refer to the MATLAB code for other parameter settings.

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.

### Version history
- Version 1.2: sven.jacobsson@ericsson.com - simplified/commented code for GitHub
- Version 1.1: studer@cornell.edu - 1-bit branch-and-bound added
- Version 1.0: sven.jacobsson@ericsson.com - minor bug fixes
- Version 0.1: sven.jacobsson@ericsson.com - simplified/commented code for GitHub
