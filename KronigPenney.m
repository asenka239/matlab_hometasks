%
%
% To launch matlab-function, you need: 
% Define k, m, a, b, U_0, Emax
%
% For example, for electron: 
%
% m = 9.109383701528 * 10^(-31);
% U_0 = -4;
% Emax = 26;
% a = 0.5*10^(-9);
% b = 2*10^(-9);
% k = [-pi/(a + b):0.01/(a + b):pi/(a + b)];
% 
% Then run: 
% E = KronigPenney(k, m, a, b, U_0, Emax);
%
% To plot, run this lines:
% plot(k,E); xlabel('k, m'); ylabel('E, eV')
% grid on

function [E] = KronigPenney(k, m, a, b, U, Emax)

% Fundamental constants:
h = 1.054571817 * 10^(-34);
eV = 1.602176634 * 10^(-19);

U_0 = U .* eV;

% Define F(E) :
function F1 = f1(E_)
    E_ = eV * E_ ;
    mu = sqrt(2*m.*E_/(h*h));
    lambda = sqrt(2*m.*(E_ - U_0)/(h*h));
    F1 = cos(mu.*a).*cos(lambda.*b)-(lambda.^2+mu.^2)./(2*mu.*lambda).*(sin(mu*a).*sin(lambda.*b));
end

% F(E) - cos(k(a+b)) to find solutions 
function  F = f(E_)
    F = f1(E_) - cos(k_*(a+b));
end

% Define small step
E = [U:10^(-3)*abs(Emax-U):Emax];
y = E;


% Now let's find intervals where our function is monotone:

for kk = 1:(length(E))
    y(kk) = f1(E(kk));
end

der = diff(y)./diff(E);
len = length(der);
extremum = U+10^(-4)*abs(Emax - U);
for ii = 1:(len - 1)
    if der(ii) * der(ii + 1) < 0
        extremum = [extremum, E(ii+1)];
    end  
end

% Great
% Now we are going to find solutions of our F(E) - cos(k(a+b)) = 0 equation:
% 
E = zeros(length(extremum), length(k));

for kk = 1:(length(extremum)-1)
    intervalE = [extremum(kk) extremum(kk + 1)];
    for ii = 1:length(k)
        k_ = k(ii);  
        E(kk,ii) = fzero(@f, intervalE);
    end
end

end
