function [x, beta, P] = MRT(s, H)
% =========================================================================
% maximal-ratio transmission (MRT) precoder
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

    % number BS antennas
    [~, B] = size(H);

    % precoding factor
    beta = sqrt(trace(H'*H))/B;
    
    % precoding matrix
    P = 1/B/beta * H';
    
    % precoded vector
    x = P*s;
                
end