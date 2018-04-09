function [x, beta] = SQUID(s,H,N0)
% =========================================================================
% squared infinity-norm relaxation with Douglas-Rachford splitting (SQUID)
%   -- inputs:
%       - s: Ux1 symbol vector
%       - H: UxB channel matrix
%       - N0: noise power spectral density (scalar)
%   -- outputs: 
%       - x: Bx1 precoded vector
%       - beta: precoding factor (scalar)
% -------------------------------------------------------------------------
% (c) 2018 Christoph Studer and Sven Jacobsson
% e-mail: studer@cornell.edu and sven.jacobsson@ericsson.com
% -------------------------------------------------------------------------
% If you use this precoder or parts of it, then you must cite our paper:
%   -- S. Jacobsson, G. Durisi, M. Coldrey, T. Goldstein, and C. Studer,
%   "Quantized precoding for massive MU-MIMO", IEEE Trans. Commun.,
%   vol. 65, no. 11, pp. 4670--4684, Nov. 2017.
% =========================================================================

    % dimensions
    [U, B] = size(H);

    % number of iterations
    iter = 50; 
    
    % gain: affects the convergence of SQUID 
    % set to 1 for large problems (e.g., 128 BS antennas) and low SNR
    % set to small values for small (ill-conditioned) problems and high SNR
    gain = 1; % default value (MUST be optimized)

    % convert to real-valued channel
    HR = [ real(H) -imag(H) ; imag(H) real(H) ];
    sR = [ real(s) ; imag(s) ];
  
    % initialize
    b = zeros(2*B,1);
    c = zeros(2*B,1);
    
    % pre-processing
    Q = HR'/((0.5/gain)*eye(2*U) + HR*HR');
    sMF = HR'*sR;
    sREG = (2*gain)*(sMF - Q*(HR*sMF));
   
    for t=1:iter % SQUID loop
        z = 2*b -c;
        a = sREG + z - Q*(HR*z);
        b = prox_infinity_norm_squared(c+a-b,2*U*B*N0);        
        c = c + a - b;
    end
    
    % extract binary solution
    x = sign(b);     
    x = 1/sqrt(2*B)*(x(1:B,1)+1i*x(B+1:2*B,1));

    % compute beta
    Hx = H*x; 
    beta = real(Hx'*s)/(norm(Hx,2)^2+U*N0);
    
    % flip negative beta
    if beta < 0
        x = -x;
        beta = -beta;
    end

end

