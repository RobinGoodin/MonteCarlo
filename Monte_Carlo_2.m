clear
clc
close
rand('seed',666);
I = 0;
for k=0:100
   I = I + 16*1/(factorial(k)*(2*k+1)^4);
end
%I = 16.2116
C = 16.1975;
f = @(x,y,z,u) exp(x*y*z*u)/sqrt(x*y*z*u);
g = @(x,y,z,u) 1 + x*y*z*u;
tic
N = 100000;
M = 0;
M2 = 0;
y = rand(N,4);
for i=1:N
    M = M + C + 16*(exp(y(i,1)^2*y(i,2)^2*y(i,3)^2*y(i,4)^2) ...
        - 1 - y(i,1)^2*y(i,2)^2*y(i,3)^2*y(i,4)^2) ;
    M2 = M2 + (16*exp(y(i,1)^2*y(i,2)^2*y(i,3)^2*y(i,4)^2))^2;
end
M = M / N;
time = toc;
M2 = M2 / (N * (N - 1));
D = M2 - M^2 / (N - 1);
S = D * time;
fprintf('N = %d\n', N);
fprintf('Origin Integral %f\n', I);
fprintf('Estimate of Integral %f\n', M);
fprintf('Elapced time %f seconds\n', time);
fprintf('Dispersion(*10^8) %f\n', D*10^8);
fprintf('Laboriousness(*10^8) %f\n', S*10^8);