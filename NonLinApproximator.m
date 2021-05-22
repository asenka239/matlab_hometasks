function [P,sgP] = NonLinApproximator(y,r,fun, P_)

% Sizes of matrix:

N = size(y,2);
M = size(P_,2);
K = size(r,1);

% Iteration parameters:

delta = 1e-6;
N_iter = 1000;


%  We will use Gauss–Newton algorithm 
% (J^T J) delta beta = J^T delta y
% So we need J and delta beta

for it=1:N_iter
    
    % Let's try to calculate Jacobian 
    
    f = zeros(N,1);
    J = zeros(N,K);
    for i=1:N
        f(i) = fun(r(:,i),P_);
        for k=1:M        
            zero = zeros(1,M);
            
            if P_(k) ~= 0        
                zero(k) = delta*abs(P_(k));
                bplus = P_ + zero; 
                bminus = P_ - zero;
                J(i,k)=0.5*(fun(r(:,i),bplus)-fun(r(:,i),bminus))/zero(k);
            else 
                zero(k) = delta;
                bplus = P_ + zero;
                J(i,k)=(fun(r(:,i),bplus)-f(i))/delta;
            end
        end
    end
    %
    % Ok, now we going to calculate delta beta:
    residial = y' - f; 
    delta_beta = (((J' * J) \ J') * residial);
    
    if sum(isnan(delta_beta)) == 0 % if everything is OK
        P_ = P_ + delta_beta';  
        if max(abs(delta_beta'./P_)) < delta
            break;
        end    
    else
        break;
    end
end

% Now, after all iterations, our P_ became true P:
P = P_;

% And the corresponding unsertanity can be again calculated as:
sgP = norm(y' - f)^2*diag(inv(J' * J))';

end
