function [x, beta] = SP(s, H, N0)
% =========================================================================
% sphere precoding (SP)
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
    
    % maximal number of iterations
    iter = 10;
    
    % dimensions
    [U, B] = size(H);
    
    % 1-bit alphabet and quantizer
    alphabet = [-1-1i; 1-1i; -1+1i; 1+1i] / sqrt(2*B);
    quantizer = @(z) (sign(real(z)) + 1i*sign(imag(z)))/sqrt(2*B);
    
    % augmented channel matrix and symbol vector
    Hb = [H; sqrt(U*N0)*eye(B)];
    sb = [s; zeros(B,1)];

    % WF precoded vector and precoding factor
    x = quantizer(WF(s, H, N0)); % initial guess
    beta = real(x'*H'*s)/(norm(H*x,2)^2+U*N0);

    % initialization
    x_list = nan(B, iter+1); x_list(:,1) = x;
    beta_list = nan(1, iter+1); beta_list(1) = beta;

    t = 1; % iteration index

        while t <= iter
            
            % QR decomposition
            if sorted_decomp == 1
                [Q,R,order] = sqr(beta_list(t)*Hb);
            else
                [Q,R] = qr(beta_list(t)*Hb); order = 1:B;
            end
            sbq = Q'*sb; 
            
            % set initial radius accoding to WF solution
            if initial_radius == 1
                Radius = norm(Q'*sb - beta_list(1)*R*x,2)^2; % radius of initial guess 
            else
                Radius = inf; % no initial radius
            end

            % initialization  
            PA = zeros(B,1); % path
            ST = zeros(B,length(alphabet)); % stack

            % add root node to stack
            Level = B; 
            ST(Level,:) = abs(sbq(Level)-R(Level,Level)*alphabet).^2;

            % begin sphere decoder
            while Level<=B   
              
                % find smallest PED in boundary    
                [minPED,idx] = min( ST(Level,:) );

                % only proceed if list is not empty
                if minPED<inf
                    
                    ST(Level,idx) = inf; % mark child as tested        
                    NewPath = [ idx ; PA(Level+1:end,1) ]; % new best path

                    % search child
                    if minPED<Radius 
                        
                        % valid candidate found
                        if  Level>1 
                        
                            % expand this best node
                            PA(Level:end,1) = NewPath;
                            Level = Level-1; % downstep
                            DF = R(Level,Level+1:end) * alphabet(PA(Level+1:end,1));
                            
                            % add to stack
                            ST(Level,:) = minPED + abs(sbq(Level)-R(Level,Level)*alphabet-DF).^2;
                        
                        else
                            
                            % valid leaf found     
                            idxhat = NewPath;
                            x = alphabet(idxhat);
                            
                            % update radius (radius reduction)
                            Radius = minPED;
                            
                        end
                        
                    end
                    
                else
                
                    % no more childs to be checked
                    Level=Level+1;
              
                end
                
            end

            % increment index counter
            t = t + 1; 

            % update precoded vector and precoding factor
            x_list(order,t) = x(:,1);
            beta_list(t) = real(x_list(:,t)'*H'*s)/(norm(H*x_list(:,t),2)^2+U*N0);

            % check stopping condition
            if sum(x_list(:,t) == x_list(:,t-1)) == B
                break;
            end

        end

        % pick solution that minimizes objective function
        [~, idxmin] = min(sum(abs(bsxfun(@minus, s, bsxfun(@times, beta_list(1:t), H*x_list(:,1:t)))).^2,1) + beta_list(1:t).^2*U*N0); 
        x = x_list(:,idxmin);
        beta = beta_list(idxmin);
        
        % flip negative beta
        if beta < 0 
            x = -x;
            beta = -beta;
        end

end


