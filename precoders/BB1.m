function [x, beta] = BB1(s, H, N0)
% =========================================================================
% 1-bit branch-and-bound (BB-1)
%   -- inputs:
%       - s: Ux1 symbol vector
%       - H: UxB channel matrix
%       - N0: noise power spectral density (scalar)
%   -- outputs:
%       - x: Bx1 precoded vector
%       - beta: precoding factor (scalar)
% -------------------------------------------------------------------------
% (c) 2018 Christoph Studer and Sven Jacobsson
% e-mail: studer@cornell.edu and sven.jacobsson@ericsson.com
% -------------------------------------------------------------------------
% If you use this precoder or parts of it, then you must cite our papers:
%   -- S. Jacobsson, G. Durisi, M. Coldrey, T. Goldstein, and C. Studer,
%   "Quantized precoding for massive MU-MIMO", IEEE Trans. Commun.,
%   vol. 65, no. 11, pp. 4670--4684, Nov. 2017.
%   -- S. Jacobsson, W. Xu, G. Durisi, and C. Studer, "MSE-optimal 1-bit 
%   precoding for multiuser MIMO via branch and bound," in Proc. IEEE Int. 
%   Conf. Acoust., Speech, Signal Process. (ICASSP), Calgary, Canada, Apr. 
%   2018, to appear.
% =========================================================================

    % tricks (true/false)
    sorted_decomp = true; % sorted QR decomposition
    initial_radius = true; % initial radius using WF solution
    prune_symmetry = true; % prune 1/2 of the tree by symmetry
    predict_future = true; % tighten bound on numerator
    
    % dimensions
    [U, B] = size(H);
    
    % 1-bit alphabet and quantizer
    alphabet = [-1-1i; 1-1i; -1+1i; 1+1i] / sqrt(2*B);
    quantizer = @(z) (sign(real(z)) + 1i*sign(imag(z)))/sqrt(2*B);

    % augmented channel matrix
    Hb = [H; sqrt(U*N0)*eye(B)];
 
    % QR decomposition
    if sorted_decomp == 1
        [~,R,order] = sqr(Hb);
    else
        [~,R] = qr(Hb); order = 1:B;
    end

    % MRT vector (scaled) before/after quantization
    zMRT = H(:,order)'*s;
    xMRT = quantizer(zMRT);

    % set initial radius accoding to WF solution
    if initial_radius == 1
        x = quantizer(WF(s,H(:,order),N0)); % initial guess (WF solution)
        radius = norm(R*x,2)^2/real(zMRT'*x)^2; % radius of initial guess 
    else
        radius = inf; % no initial radius
    end

    % initialization
    PA = zeros(B,1); % path
    ST = zeros(B,length(alphabet)); % stack
    
    % preprocessing
    num_next = nan(4,B); % numerator (next step)
    num_potential = cell(1,B); % numerator (potential steps)
    den_present = nan(4,B); % denominator (present)
    den_future = nan(1,B); % denominator (future)
    for l = 1:B
        num_potential{:,l} = R(1:l-1,l)*alphabet.';
        num_next(:,l) = R(l,l)*alphabet;
        den_present(:,l) = real(zMRT(l)'*alphabet);
        den_future(l) = real(zMRT(1:l-1)'*xMRT(1:l-1));
    end
    
    % predicting the future
    if predict_future == 1
        iRt = cell(1,B);
        eigmin = nan(1,B);
        for l=2:B
            iRt{l} = inv(R(1:l-1,1:l-1));
            eigmin(l) = min(svd(R(1:l-1,1:l-1)))^2;
        end
    end
    
    % root node
    level = B;
    
    % numerator (lower bound)
    if predict_future == 1
        b = iRt{level}*(R(1:level-1,level)*alphabet.');        
        num_future = eigmin(level)*sum(abs(quantizer(b)-b).^2,1).';        
        numerator = abs(num_next(:,l)).^2 + num_future;
    else
        numerator = abs(num_next(:,l)).^2; % present only
    end
    
    % denominator (upper bound)
    denominator = (abs(den_present(:,level)) + abs(den_future(level))).^2; % present and future
    
    % add root node to stack
    ST(level,:) = numerator./denominator;
    
    % prune half of the tree by symmetry
    if prune_symmetry == 1
        ST(level,3:4) = inf; % exclude symmetric solutions
    end
    
    % sphere decoding
    while level <= B
        
        % find smallest PED in boundary
        [minPED,idx] = min(ST(level,:));
        
        % proceed only if list is not empty
        if minPED<inf
            
            ST(level,idx) = inf; % mark child as tested
            NewPath = [idx; PA(level+1:end,1)]; % new best path
            
            % search child
            if minPED < radius
               
                % valid candidate found
                if level>1 
                    
                    % expand this best node
                    PA(level:end,1) = NewPath;
                    level = level - 1; % downstep
                    
                    % numerator (lower bound)
                    num_past = norm(R(level+1:end,level+1:end)*alphabet(PA(level+1:end,1)),2)^2;
                    num_present = abs(num_next(:, level) + R(level,level+1:end) * alphabet(PA(level+1:end,1))).^2;
                    if predict_future == 1 && level>1
                        b = iRt{level}*(num_potential{:,level} + (R(1:level-1,level+1:end)*alphabet(PA(level+1:end,1)))*ones(1,length(alphabet)));
                        num_future = eigmin(level)*sum(abs(b-quantizer(b)).^2,1).';
                    else
                       num_future = 0; 
                    end
                    numerator =  num_past + num_present + num_future;
                    
                    % denominator (upper bound)
                    den_past = real(zMRT(level+1:end,1)'*alphabet(PA(level+1:end,1)));
                    denominator = ( abs(den_present(:,level) +  den_past) +  abs(den_future(level)) ).^2;
                    
                    % add to stack
                    ST(level,:) = numerator./denominator;
                    
                else
                    
                    % valid leaf found
                    idxhat = NewPath;
                    x = alphabet(idxhat);
                    
                    % update radius (radius reduction)
                    radius = minPED;
                    
                end
                
            end
            
        else
            
            % no more child nodes to be checked
            level=level+1;
            
        end
        
    end
    
    % re-order vector and compute beta
    x(order,1) = x;
    beta = real(x'*H'*s)/(norm(H*x,2)^2+U*N0);
    
    % flip negative beta
    if beta < 0 
        x = -x;
        beta = -beta;
    end
    
end

