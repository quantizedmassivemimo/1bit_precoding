function xk = prox_infinity_norm_squared(w,lambda)
% proximal operator for the squared infinity norm.

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