function Q = ElectroStaticBalls(XYZ,R,F)

N = size(F);

% First we will sheck that our balls do not overlap:

for i = 1:N
    for j = (i + 1):N
        if norm(XYZ(:,i)-XYZ(:,j)) <= R(i) + R(j)
            error('Overlap!')
        end
    end
end

% OK, now we can get our matrix to solve System of linear equations:

Matrix = zeros(N);

for i = 1:N
    % For the diagonal elements:
    Matrix(i, i) = 1 / R(i);
    
    % Solve: 
    for j = 1:N
        if(i ~= j)
            Matrix(i, j) = 1 / norm(XYZ(:,i)-XYZ(:,j));
        end
    end
end

% WIth knowledge of potential F we obtain our charges Q:
Q = Matrix \ F;