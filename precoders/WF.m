function [x, beta, P] = WF(s, H, N0)
% =========================================================================
% Wiener filter (WF) precoder
%   -- inputs:
%       - s: Ux1 complex-valued symbol vector
%       - H: UxB complex-valued channel matrix
%       - N0: noise power spectral density (scalar)
%   -- outputs: 
%       - x: Bx1 complex-valued precoded vector
%       - beta: precoding factor (scalar)
%       - P: BxU precoding matrix
% -------------------------------------------------------------------------
% (c) 2017 Christoph Studer and Sven Jacobsson
% e-mail: studer@cornell.edu and sven.jacobsson@ericsson.com
% =========================================================================

    % number of UEs
    [U, ~] = size(H);
    
    % precoding matrix (before normalization)
    T = H' / (H*H' + U*N0*eye(U));

    % precoding factor
    beta = sqrt(real(trace((T*T'))));
    
    % precoding matrix
    P = 1/beta*T;
    
    % precoded vector
    x = P*s;

end
