% =========================================================================
% Sphere precoding (SP)
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

function [x, beta] = SP(par,s,H,N0)

    MAX = 10; % maximal number of iterations
    
    if par.L == 2

        % -- preprocessing

        % initialization
        beta_list = nan(1, MAX+1);
        x_list = nan(par.B, MAX+1);

        % append symbol vector and channel matrix
        Hb = [H; sqrt(par.U*N0)*eye(par.B)];
        sb = [s; zeros(par.B,1)];

        % ZF precoded vector and precoding factor
        x_list(:,1) = par.quantizer(ZF(s, H));
        beta_list(1) = real(x_list(:,1)'*H'*s)/(norm(H*x_list(:,1),2)^2+par.U*N0);


        t = 1; % iteration index

        while t <= MAX

          % -- initialization  
          Radius = inf;
          PA = zeros(par.B,1); % path
          ST = zeros(par.B,length(par.alphabet)); % stack  

          [Q,R] = qr(beta_list(t)*Hb,0);  
          shat = Q'*sb;    

          % -- add root node to stack
          Level = par.B; 
          ST(Level,:) = abs(shat(Level)-R(Level,Level)*par.alphabet.').^2;

          % -- begin sphere decoder
          while ( Level<=par.B )          
            % -- find smallest PED in boundary    
            [minPED,idx] = min( ST(Level,:) );

            % -- only proceed if list is not empty
            if minPED<inf
              ST(Level,idx) = inf; % mark child as tested        
              NewPath = [ idx ; PA(Level+1:end,1) ]; % new best path

              % -- search child
              if ( minPED<Radius )
                % -- valid candidate found
                if ( Level>1 )                  
                  % -- expand this best node
                  PA(Level:end,1) = NewPath;
                  Level = Level-1; % downstep
                  DF = R(Level,Level+1:end) * par.alphabet(PA(Level+1:end,1)).';
                  ST(Level,:) = minPED + abs(shat(Level)-R(Level,Level)*par.alphabet.'-DF).^2;
                else
                  % -- valid leaf found     
                  idxhat = NewPath;
                  x_temp = par.alphabet(idxhat).';
                  % -- update radius (radius reduction)
                  Radius = minPED;    
                end
              end      
            else
              % -- no more childs to be checked
              Level=Level+1;      
            end    
          end

          t = t + 1; % increment index counter;

          % update precoded vector and 
          x_list(:,t) = x_temp;
          beta_list(t) = real(x_temp'*H'*s)/(norm(H*x_temp,2)^2+par.U*N0);

          if sum(x_list(:,t) == x_list(:,t-1)) == par.B
              break;
          end

        end

        % -- pick solution that minimizes objective function
        J = sum(abs(bsxfun(@minus, s, bsxfun(@times, beta_list(1:t), H*x_list(:,1:t)))).^2,1) + beta_list(1:t).^2*par.U*N0;
        [~, idxmin] = min(J); 

        x = x_list(:,idxmin);
        beta = beta_list(idxmin);
        
    else
        
        error('SP: only 1-bit DACs supported!');
        
    end
    
    if beta < 0
        warning('!!!!');
    end
  
end