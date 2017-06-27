% =========================================================================
% Squared infinity-norm relaxation with Douglas-Rachford splitting (SQUID)
%   -- inputs:
%       - par: struct of simulation parameters
%       - s: Ux1 complex-valued symbol vector
%       - H: UxB complex-valued channel matrix
%       - N0: noise power spectral density (scalar)
%   -- outputs: 
%       - x: Bx1 complex-valued precoded vector
%       - beta: precoding factor (scalar)
% -------------------------------------------------------------------------
% (c) 2017 Christoph Studer and Sven Jacobsson
% e-mail: studer@cornell.edu and sven.jacobsson@ericsson.com
% =========================================================================

function [x, beta] = SQUID(par,s,H,N0)
    
    % convert to real-valued channel
    HR = [ real(H) -imag(H) ; imag(H) real(H) ];
    sR = [ real(s) ; imag(s) ];
  
    % initialize
    x = zeros(par.B*2,1);
    y = zeros(par.B*2,1);
    
    gain = 1; % set to 1 for large problems; small values for small, ill-conditioned problems
    epsilon = 1e-5; % needs to be about 1e-5 for accurate results
    
    % pre-processing
    A = HR'*HR + 0.5/gain*eye(par.B*2);
    sREG = A\(HR'*sR);
    
    if par.L == 2
      
        % SQUID loop
        for t=1:100
            u = sREG + 0.5/gain*(A\(2*x-y));
            xold = x;
            x = prox_infinityNorm2(y+u-x,2*2*par.U*par.B*N0);        
            if norm(x-xold)/norm(x)<epsilon
                break;
            end
            y = y + u - x;
        end
        
    else
        
        error('SQUID: only 1-bit DACs supported!');
        
    end
    
    % extract binary solution
    xRest = sign(x);        
    x = 1/sqrt(2*par.B)*(xRest(1:par.B,1)+1i*xRest(par.B+1:2*par.B,1));
    
    % compute output gains
    beta = real(x'*H'*s)/(norm(H*x,2)^2+par.U*N0);
    
    % check (and fix) if beta is negative
    if beta < 0
        x = -x;
        beta = -beta;
    end

end

% proximal mapping of the infinity-norm-squared.
% perform prox operator: min lambda*||x||_inf^2 + ||x-w||^2
function xk = prox_infinityNorm2(w,lambda)

    N = length(w);
    wabs = abs(w);
    ws = (cumsum(sort(wabs,'descend')))./(lambda+(1:N)');
    alphaopt = max(ws);
    
    if alphaopt>0 
      xk = min(wabs,alphaopt).*sign(w); % truncation step
    else
      xk = zeros(size(w)); % if t is big, then solution is zero
    end   
    
end

