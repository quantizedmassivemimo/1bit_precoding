function [Qa,R,order] = sqr(H)
% sorted QR Decomposition

    % initializaion
    Ntx = size(H,2);
    Nrx = size(H,1);
    Q=H;
    R=zeros(Ntx);
    P=eye(Ntx);
    
    % sorting
    for i= 1:Ntx
        
        s=diag(Q(:,i:Ntx)'*Q(:,i:Ntx));    
        [~,idx]=min(s); 
        idx=idx+i-1;
        tmp=Q(:,i); Q(:,i)=Q(:,idx); Q(:,idx)=tmp;
        tmp=P(:,i); P(:,i)=P(:,idx); P(:,idx)=tmp;
        tmp=R(:,i); R(:,i)=R(:,idx); R(:,idx)=tmp;
        R(i,i)=sqrt(s(idx-i+1));
        Q(:,i)=Q(:,i)/R(i,i);
        
        for k= i+1:Ntx
            
            R(i,k)=Q(:,i)'*Q(:,k);
            Q(:,k)=Q(:,k)-R(i,k)*Q(:,i);
            
        end
        
    end
    
    Qa    = Q(1:Nrx,1:Ntx); % take upper left part only
    R     = R(1:Ntx,:); % take upper half only (lower half is zero)
    order = (1:Ntx)*P; % order of sorted decomposition
    
end

