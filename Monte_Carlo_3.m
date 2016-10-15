clear
clc
close
rand('seed',666);
I = 0;
for k=0:10
   I = I + 16*1/(factorial(k)*(2*k+1)^4);
end
%I = 16.2116
A = 81/1312;
C = 16.1975;
f = @(x,y,z,u) A/C*(1+x*y*z*u);
g = @(x,y,z,u) 1 + x*y*z*u;
tic
N = 100000;
M = 0;
M2 = 0;
p1 = 2^4*A;
p2 = (2/3)^4*A;
g = rand(N,5);
for i=1:N
    if( g(i,5) < p1)
        x = g(i,1)^2;
        y = g(i,2)^2;
        z = g(i,3)^2;
        u = g(i,4)^2;
    else
        x = g(i,1)^(2/3);
        y = g(i,2)^(2/3);
        z = g(i,3)^(2/3);
        u = g(i,4)^(2/3);
    end
    M = M + exp(x^2*y^2*z^2*u^2) / (A*(1+x*y*z*u));
    M2 = M2 + (exp(x^2*y^2*z^2*u^2) / (A*(1+x*y*z*u)))^2;   
end
M = M / N;
time = toc;
M2 = M2 / (N * (N - 1));
D = M2 - M^2 / (N - 1);
S = D * time;
fprintf('Origin Integral %f\n', I);
fprintf('Estimate of Integral %f\n', M);
fprintf('Elapced time %f seconds\n', time);
fprintf('Dispersion(*10^8) %f\n', D*10^8);
fprintf('Laboriousness(*10^8) %f\n', S*10^8);