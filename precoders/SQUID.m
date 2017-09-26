function [x, beta] = SQUID(par,s,H,N0)
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
    
    % convert to real-valued channel
    HR = [ real(H) -imag(H) ; imag(H) real(H) ];
    sR = [ real(s) ; imag(s) ];
  
    % initialize
    b = zeros(2*par.B,1);
    c = zeros(2*par.B,1);
    
    gain = 1; % set to 1 for large problems; small values for small, ill-conditioned problems
    iter = 50; % number of iterations
    
    % pre-processing
    Q = HR'/((0.5/gain)*eye(2*par.U) + HR*HR');
    sMF = HR'*sR;
    sREG = (2*gain)*(sMF - Q*(HR*sMF));
   
    if par.L == 2
        for t=1:iter % SQUID loop
            z = 2*b -c;
            a = sREG + z - Q*(HR*z);
            b = prox_infinityNorm2(c+a-b,2*par.U*par.B*N0);        
            c = c + a - b;
        end
    else
        error('SQUID: only 1-bit DACs supported!');
    end
    
    % extract binary solution
    x = sign(b);     
    x = 1/sqrt(2*par.B)*(x(1:par.B,1)+1i*x(par.B+1:2*par.B,1));

    % compute beta
    Hx = H*x; 
    beta = real(Hx'*s)/(norm(Hx,2)^2+par.U*N0);
    
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
    ws = (cumsum(sort(wabs,'descend')))./(2*lambda+(1:N)');
    alphaopt = max(ws);
    
    if alphaopt>0 
      xk = min(wabs,alphaopt).*sign(w); % truncation step
    else
      xk = zeros(size(w)); % if t is big, then solution is zero
    end   
    
end

