function [x_sdr1, beta_sdr1, x_sdrr, beta_sdrr] = SDR(par,s,H,N0)
% =========================================================================
% Semidefinite relaxation (SDR)
%   -- inputs:
%       - par: struct of simulation parameters
%       - s: Ux1 complex-valued symbol vector
%       - H: UxB complex-valued channel matrix
%       - N0: noise power spectral density (scalar)
%   -- outputs: 
%       - x_sdr1: Bx1 complex-valued quantized precoded vector (rank-one approx.)
%       - beta_sdr1: precoding factor (rank-one approx.)
%       - x_sdrr: Bx1 complex-valued quantized precoded vector (randomization)
%       - beta_sdrr: precoding factor (randomization)
% -------------------------------------------------------------------------
% (c) 2017 Christoph Studer and Sven Jacobsson
% e-mail: studer@cornell.edu and sven.jacobsson@ericsson.com
% =========================================================================

    MAX = 1e2; % number of candidate solutions (randomizations)

    % convert to real-valued channel
    sR = [real(s); imag(s)];
    HR = [real(H), -imag(H); imag(H), real(H)];
    
    if par.L == 2
        
        % solve SDP
        cvx_begin quiet
            TR = [HR'*HR+par.U*N0*eye(2*par.B), -HR'*sR; -sR'*HR, norm(sR,2)^2];
            variable BR(2*par.B+1,2*par.B+1) symmetric
            minimize trace(TR*BR)
            subject to
            for j = 2:2*par.B 
                BR(1,1) == BR(j,j);
            end
            BR(2*par.B+1,2*par.B+1) == 1;
            BR == semidefinite(2*par.B+1);
        cvx_end
        
    else
        
        error('SDR: only 1-bit DACs supported!');
        
    end
    
    % eigenvalue decomposition
    [V, D] = eig(BR);
    
    % -- Rank-one approximation
    
    % find maximum eigenvalue 
    [~, idxmax] = max(diag(D)); 

    % quantize to feasible solution and convert from real ro complex
    xR_sdr1 = par.quantizer(sqrt(D(idxmax,idxmax))*V(:,idxmax));
    x_sdr1 = sign(xR_sdr1(end)) * (xR_sdr1(1:par.B,1)+1i*xR_sdr1(par.B+1:2*par.B,1));

    % compute precoding factor
    beta_sdr1 = real(x_sdr1'*H'*s)/(norm(H*x_sdr1,2)^2+par.U*N0);

    % -- Randomization
    
    D = diag(max(diag(D), 0)); % check eigenvalues
    BR = V*D*V'; % ensures that matrix is positive semidefinite
    
    % generate candidate vectors with correspodnding precoding factor
    xR = [xR_sdr1, par.quantizer(mvnrnd(zeros(2*par.B+1,1), BR , MAX)')]; 
    HxR = HR*xR(1:end-1,:);
    beta = sR'*HxR./(dot(HxR,HxR) + par.U*N0);

    % compute obj. func. for all candidates 
    J = sum(abs(bsxfun(@minus, sR, bsxfun(@times, beta, HxR))).^2,1) + beta.^2*par.U*N0;
 
    % pick candidate that maximizes obj. function
    [~, idxmax] = min(J);

    % real ro complex
    xR_sdrr = xR(:,idxmax); 
    x_sdrr = sign(xR_sdrr(end)) * (xR_sdrr(1:par.B,1)+1i*xR_sdrr(par.B+1:2*par.B,1));

    % compute precoding factor
    beta_sdrr = real(x_sdr1'*H'*s)/(norm(H*x_sdr1,2)^2+par.U*N0); 
    
end

