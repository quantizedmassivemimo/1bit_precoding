% =========================================================================
% Exhaustive search (EXS) precoder
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

function [x, beta] = EXS(par,s,H,N0)

    par.EXS.type = 'complex'; % real or complex 

    if par.L == 2

        % all possible quantization outcomes
        x_list = par.alphabet;
        for i = 1:par.B-1
            x_list = combvec(x_list, par.alphabet);
        end

        % consider only outcomes with positive gain
        x_list = x_list(:, real(s'*H*x_list)>0); 

        if strcmpi(par.EXS.type, 'complex')

            % received noiseless signal
            Hx_list = H*x_list; 

            % precoding factor
            beta = real(s'*Hx_list)./(dot(Hx_list,Hx_list) + par.U*N0); 

            % objective function
            J = sum(abs(bsxfun(@minus, s, bsxfun(@times, beta, Hx_list))).^2,1) + beta.^2*par.U*N0;

        elseif strcmpi(par.EXS.type, 'real')

            % convert to real-valued channel
            HR = [real(H), -imag(H); imag(H), real(H)];
            sR = [real(s); imag(s)]; 
            xR_list = [real(x_list); imag(x_list)];

            % received noiseless signal
            HxR_list = HR*xR_list; 

            % precoding factor
            beta = sR'*HxR_list./(dot(HxR_list,HxR_list) + par.U*N0); 

            % objective function
            J = sum(abs(bsxfun(@minus, sR, bsxfun(@times, beta, HxR_list))).^2,1) + beta.^2*par.U*N0;

        end

        % pick solution that minimizes objective function
        [~, idxmin] = min(J); 
        x = x_list(:,idxmin);
        beta = beta(idxmin);
        
    else
        
        error('EXS: only 1-bit DACs supported!');
        
    end
    
end

