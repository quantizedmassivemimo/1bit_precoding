function [x, beta, P] = ZF(s, H)
% =========================================================================
% Zero-forcing (ZF) precoder
%   -- inputs:
%       - s: Ux1 complex-valued symbol vector
%       - H: UxB complex-valued channel matrix
%   -- outputs: 
%       - x: Bx1 complex-valued precoded vector
%       - beta: precoding factor (scalar)
%       - P: BxUprecodin matrix
% -------------------------------------------------------------------------
% (c) 2017 Christoph Studer and Sven Jacobsson
% e-mail: studer@cornell.edu and sven.jacobsson@ericsson.com
% =========================================================================

    % precoding factor
    beta = sqrt(trace((H*H')^-1)); 
    
    % precoding matrix
    P = 1/beta * H'/(H*H');
    
    % precoded vector
    x = P*s;

end
