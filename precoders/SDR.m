function [x_sdr1, beta_sdr1, x_sdrr, beta_sdrr] = SDR(s,H,N0)
% =========================================================================
% semidefinite relaxation (SDR)
%   -- inputs:
%       - s: Ux1 symbol vector
%       - H: UxB channel matrix
%       - N0: noise power spectral density (scalar)
%   -- outputs: 
%       - x_sdr1: Bx1 precoded vector (rank-one approx.)
%       - beta_sdr1: precoding factor (rank-one approx.)
%       - x_sdrr: Bx1 precoded vector (randomization)
%       - beta_sdrr: precoding factor (randomization)
% -------------------------------------------------------------------------
% (c) 2018 Christoph Studer and Sven Jacobsson
% e-mail: studer@cornell.edu and sven.jacobsson@ericsson.com
% -------------------------------------------------------------------------
% If you use this precoder or parts of it, then you must cite our paper:
%   -- S. Jacobsson, G. Durisi, M. Coldrey, T. Goldstein, and C. Studer,
%   "Quantized precoding for massive MU-MIMO", IEEE Trans. Commun.,
%   vol. 65, no. 11, pp. 4670--4684, Nov. 2017.
% =========================================================================

    MAX = 1e2; % number of candidate solutions (randomizations)
    
    % dimensions
    [U, B] = size(H);
    
    % 1-bit quantizer
    quantizer = @(z) (sign(real(z)) + 1i*sign(imag(z)))/sqrt(2*B);

    % convert to real-valued channel
    sR = [real(s); imag(s)];
    HR = [real(H), -imag(H); imag(H), real(H)];
    
    % solve SDP
    cvx_begin quiet
        TR = [HR'*HR+U*N0*eye(2*B), -HR'*sR; -sR'*HR, norm(sR,2)^2]; %#ok<NASGU>
        variable BR(2*B+1,2*B+1) symmetric
        minimize trace(TR*BR)
        subject to
        for j = 2:2*B 
            BR(1,1) == BR(j,j); %#ok<NODEF,EQEFF>
        end
        BR(2*B+1,2*B+1) == 1; %#ok<EQEFF>
        BR == semidefinite(2*B+1); %#ok<EQEFF>
    cvx_end
        
    % eigenvalue decomposition
    [V, D] = eig(BR); 
    
    % -- rank-one approximation
    
    % find maximum eigenvalue 
    [~, idxmax] = max(diag(D)); 

    % quantize to feasible solution and convert from real ro complex
    xR_sdr1 = quantizer(sqrt(D(idxmax,idxmax))*V(:,idxmax));
    x_sdr1 = sign(xR_sdr1(end)) * (xR_sdr1(1:B,1)+1i*xR_sdr1(B+1:2*B,1));

    % compute precoding factor
    beta_sdr1 = real(x_sdr1'*H'*s)/(norm(H*x_sdr1,2)^2+U*N0);

    % -- randomization
    
    D = diag(max(diag(D), 0)); % check eigenvalues
    BR = V*D*V'; % ensures that matrix is positive semidefinite
    
    % generate candidate vectors with correspodnding precoding factor
    xR = [xR_sdr1, quantizer(mvnrnd(zeros(2*B+1,1), BR , MAX)')]; 
    HxR = HR*xR(1:end-1,:);
    beta = sR'*HxR./(dot(HxR,HxR) + U*N0);

    % compute objective function for all candidates 
    J = sum(abs(bsxfun(@minus, sR, bsxfun(@times, beta, HxR))).^2,1) + beta.^2*U*N0;
 
    % pick candidate that maximizes obj. function
    [~, idxmax] = min(J);

    % real to complex
    xR_sdrr = xR(:,idxmax); 
    x_sdrr = sign(xR_sdrr(end)) * (xR_sdrr(1:B,1)+1i*xR_sdrr(B+1:2*B,1));

    % compute precoding factor
    beta_sdrr = real(x_sdr1'*H'*s)/(norm(H*x_sdr1,2)^2+U*N0); 
    
end

