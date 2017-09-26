function [x, beta, P] = MRT(s, H)
% =========================================================================
% Maximal-ratio transmission (MRT) precoder
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

    % number BS antennas
    [~, B] = size(H);

    % precoding factor
    beta = sqrt(trace(H'*H))/B;
    
    % precoding matrix
    P = 1/B/beta * H';
    
    % precoded vector
    x = P*s;
                
end