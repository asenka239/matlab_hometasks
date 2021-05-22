function [P,sgP] = LinApproximator(y,r,funcs)

% Sizes of matrix:

M = size(funcs,2);
N = size(y,2);
K = size(r,1);

% Let's now define phi = funcs(r)

phi = zeros(N,M);
for i=1:N
    for j=1:M
        phi(i,j) = funcs{j}(r(:,i));
    end
end

% Linear least squares methods can be solved (in terms of linear system) 
% P = (phi^T * phi) ^-1 * phi^T * y

P = inv(phi' * phi)  * phi' * y';

% And the corresponding unsertanity can be easily calculated as:

sgP = norm(y' - phi*P)^2 * diag(inv(phi' * phi))';

end
