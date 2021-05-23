function [F,X,Y,P] = SpherePotential(XYZ,Q,R,r0,a,b,Dx,Dy,Nxy)
N = length(Q);

% First we need to make up the coordinates
% We need to get projection of b to y axis:
b = b - a .* dot(a,b) ./ dot(a,a);
% And our Transition Matrix:
P = [a(1) b(1); 
     a(2) b(2); 
     a(3) b(3)];

% Now we can make up with the nodes: 
X = [Dx(1) : (Dx(2) - Dx(1)) / (Nxy(2) - 1) : Dx(2)];
Y = [Dy(1) : (Dy(2) - Dy(1)) / (Nxy(1) - 1) : Dy(2)];

% And here we go to potential calculation:
F = zeros(Nxy(1),Nxy(2)); % potentials

for i = 1:Nxy(2)     % loop over
    for j = 1:Nxy(1) % all nodes
        
        % This is our r0 in new coordinates
        r0_ = [r0(1) + a(1) * X(i) + b(1) * Y(j); 
               r0(2) + a(2) * X(i) + b(2) * Y(j); 
               r0(3) + a(3) * X(i) + b(3) * Y(j)];
           
        for k = 1:N % loop over charges 
            
            if(norm(XYZ(:,k)-r0_) > R(k)) % node not in ball
                F(j,i) = F(j, i) + Q(k) / norm(XYZ(:,k)-r0_);

            else                          % node in ball
                F(j,i) = F(j, i) + Q(k) / R(k); 
            end
        end
    end
end      

% Finally, making grid:

X = repmat(X, Nxy(1), 1);
Y = repmat(Y', 1, Nxy(2));

end

