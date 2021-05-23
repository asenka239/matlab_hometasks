function [Q,D] = ElectroStaticDipoles(XYZ,R,F)
XYZ = XYZ';
N = length(F);
% Adding the electric fields 3N-components: 
F = [F; zeros(3 * N, 1)];

% First we will sheck that our balls do not overlap:

for i_ = 1:N
    for j_ = (i_ + 1):N
        if norm(XYZ(:,i_)-XYZ(:,j_)) <= R(i_) + R(j_)
            error('Overlap!')
        end
    end
end

% OK, now we can get our matrix to solve System of linear equations:
Matrix = zeros(4 * N);
%
% 
% We need some conditions and calculations for our potential: 
%
for i = 1 : N
    for j = 1 : 4 * N
        j_ = mod(j, N);
        
        if(j_ == 0)
            j_ = N;
        end
        
        coordinates = (j - j_) / N;
        r_ = XYZ(1:3,i) - XYZ(1:3,j_);
        
        if(i == j)
            Matrix(i, j) = 1 / R(i);
        else
            if(coordinates == 0)
                Matrix(i, j) = 1 / norm(r_);
            else
                if(j_ == i)
                    Matrix(i, j) = 0;
                else
                    Matrix(i, j) = r_(coordinates) / norm(r_) ^ 3;
                end
            end
        end   
    end
end

%
% End now we need to make up with the E_field inside the balls:

for i = (N + 1) : (4 * N)
    for j = 1 : (4 * N)
        j_ = mod(j, N);
        i_ = mod(i, N);
        
        if(j_ == 0)
            j_ = N;
        end
        
        if(i_ == 0)
            i_ = N;
        end
        
        coord_p = (j - j_) / N;
        coord_E = (i - i_) / N;
        r_ = XYZ(1:3,i_) - XYZ(1:3,j_);
        if(i_ == j)
            Matrix(i, j) = 0;
        else
            if(i == j)
                Matrix(i, j) = 1 / R(i_) ^ 3;
            else
                if(coord_p == 0)
                    Matrix(i, j) = r_(coord_E) / norm(r_) ^ 3;
                else
                    if(j_ == i_)
                        Matrix(i, j) = 0;
                    else
                        if(coord_E == coord_p)
                            Matrix(i, j) = (3 * r_(coord_E) ^ 2 - norm(r_)^2 ) / norm(r_) ^ 5;
                        else
                            Matrix(i, j) = 3 * r_(coord_E) * r_(coord_p)/ norm(r_) ^ 5;
                        end
                    end
                end
            end
        end 
    end
end

% With knowledge of potential and fileds F 
% we obtain our charges Q and dipole D
% We need temp var and then use its projection:
temp_Q = Matrix \ F;

Q = temp_Q(1:N, 1);
D = -[temp_Q(N+1:2*N,1)'; 
      temp_Q(2*N+1:3*N,1)'; 
      temp_Q(3*N+1:4*N,1)']';

end



