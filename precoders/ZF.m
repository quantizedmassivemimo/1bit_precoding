function [x, beta, P] = ZF(s, H)
% =========================================================================
% zero-forcing (ZF) precoder
%   -- inputs:
%       - s: Ux1 symbol vector
%       - H: UxB channel matrix
%   -- outputs: 
%       - x: Bx1 precoded vector
%       - beta: precoding factor (scalar)
%       - P: BxU precoding matrix
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
