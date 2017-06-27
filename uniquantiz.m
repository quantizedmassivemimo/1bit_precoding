% =========================================================================
% Uniform quantizer
%   -- inputs (optional):
%       - y: Mx1 matrix of complex-valued quantizer input
%       - lsb: least-significant bit (step size)
%       - L: number of quantization levels
%   -- outputs: 
%       - v: Mx1 matrix of complex-valued quantizer output
%       - e: Mx1 matrix of complex-valued quantization error
%       - vl: vector with L quantization labels
%       - vt: vector with L+1 quantization thresholds
%       - clip: clip level of the quantizer
% -------------------------------------------------------------------------
% (c) 2017 Christoph Studer and Sven Jacobsson
% e-mail: studer@cornell.edu and sven.jacobsson@ericsson.com
% =========================================================================


function [v, e, vl, vt, clip] = uniquantiz(y, lsb, L)

    % default parameters (to illustrate quantizer input-output relation)
    if nargin == 0
        L = 7;
        lsb = 2;
        y = .75*(linspace(-L*lsb, L*lsb, 1e4) + 1i*linspace(-L*lsb, L*lsb, 1e4));
    end
    
    % set clip level
    clip = lsb*L/2;
    
    % clip signal
    if isreal(y)
        yc = max(min(y,clip-lsb/1e5),-(clip-lsb/1e5));
    else
        yc = max(min(real(y),clip-lsb/1e5),-(clip-lsb/1e5)) + 1i*max(min(imag(y),clip-lsb/1e5),-(clip-lsb/1e5));
    end
    
    % quantizer
    if mod(L,2) == 0
        Q = @(x) lsb*floor(x/lsb) + lsb/2; % midrise quantizer (without clipping)
    else
        Q = @(x) lsb*floor(x/lsb + 1/2); % midtread quantizer (without clipping)
    end
    
    % quantize signal
    if isreal(y)
        v = Q(yc);
    else
        v = Q(real(yc)) + 1i*Q(imag(yc));
    end
    
    % quantization error
    e = v - y;
    
    % uniform quantization labels
    vl = lsb *((0:L-1) - (L-1)/2);
    
    % uniform quantization thresholds
    vt = [-10^100, bsxfun(@minus, vl(:,2:end), lsb/2), 10^100];	
    
    if nargin == 0
        
        figure(1); clf; % illustrate quantizer
        
        maxamp = max(max(real(y), real(v)));
        
        % show quantizater-mapping function (real-valued input)
        subplot(1,2,1); hold all;
        plot(real(y), real(y),'k--');
        stairs(real(y), real(v),'r-');
        line([-maxamp,maxamp],[0,0],'color','k');
        line([0,0],[-maxamp,maxamp],'color','k');
        axis equal;box on;
        axis([-maxamp,maxamp,-maxamp,maxamp]);
        xlabel('input signal','fontsize',14);
        ylabel('quantized signal','fontsize',14);
        set(legend('input', 'quant.'), 'location', 'northwest','fontsize',12);
        
        % show quantization error (real part)
        subplot(1,2,2); hold all;
        plot(real(y), real(y),'k--');
        plot(real(y), real(e),'r-');
        line([-maxamp,maxamp],[0,0],'color','k');
        line([0,0],[-maxamp,maxamp],'color','k');
        axis equal;box on;
        axis([-maxamp,maxamp,-maxamp,maxamp]);
        xlabel('input signal','fontsize',14);
        ylabel('quantization error','fontsize',14);
        set(legend('input', 'quant. error'), 'location', 'northwest','fontsize',12);
        
  
    end

end








