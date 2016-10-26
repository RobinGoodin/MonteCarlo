clc
clear
close
%% Задание основных параметров
set(0,'RecursionLimit',10000)
rand('seed',666);
N = 100;
N_start  = 10;
N_step = 10;
func = @(x,y,z,u) exp(x*y*z*u)/sqrt(x*y*z*u);
f = @(x,y,z,u) 16*exp(x*y*z*u);
p = @(x,y,z,u) 1/(16*sqrt(x*y*z*u));
I = 0;
for k=0:1000
   I = I + 16*1/(factorial(k)*(2*k+1)^4);
end

%% Вычисление последовательности Холтона
M = N;
p2 = zeros(1, N);
p3 = zeros(1, N);
p5 = zeros(1, N);
for k = 1:N
    a2 = zeros(1, M);
    a3 = zeros(1, M);
    a5 = zeros(1, M);
    
    a2(1) = mod(k, 2);
    p2(k) = a2(1)/2;
    
    a3(1) = mod(k, 3);
    p3(k) = a3(1)/3;

    a5(1) = mod(k, 5);
    p5(k) = a5(1)/5;
    
    for j = 2:k+5
        a2(j) = mod(rec_fix(j-1, k, 2), 2);
        p2(k) = p2(k) + a2(j)*2^(-j);
        
        a3(j) = mod(rec_fix(j-1, k, 3), 3);
        p3(k) = p3(k) + a3(j)*3^(-j);
        
        a5(j) = mod(rec_fix(j-1, k, 5), 5);
        p5(k) = p5(k) + a5(j)*5^(-j);
    end
    if(mod(k,100) == 0)
        k
    end
end

%% Старый алгоритм
tic
M = 0;
% M2 = 0;
y = rand(N,4);
for i=1:N
    M = M + 16*exp(y(i,1)^2*y(i,2)^2*y(i,3)^2*y(i,4)^2);
    M2 = M2 + (16*exp(y(i,1)^2*y(i,2)^2*y(i,3)^2*y(i,4)^2))^2;
end
M = M / N;
% time = toc;
% M2 = M2 / (N * (N - 1));
% D = M2 - M^2 / (N - 1);
% S = D * time;
% fprintf('N = %d\n', N);
% fprintf('Origin Integral %f\n', I);
% fprintf('Estimate of Integral %f\n', M);
% fprintf('Elapced time %f seconds\n', time);
% fprintf('Dispersion(*10^8) %f\n', D*10^8);
% fprintf('Laboriousness(*10^8) %f\n', S*10^8);

%% Новый алгоритм
% tic
M = 0;
% M2 = 0;
for i=1:N
    M = M + 16*exp((i/N)^2*p2(i)^2*p3(i)^2*p5(i)^2);
    M2 = M2 + (16*exp((i/N)^2*p2(i)^2*p3(i)^2*p5(i)^2))^2;
end
M = M / N;
% time = toc;
% M2 = M2 / (N * (N - 1));
% D = M2 - M^2 / (N - 1);
% S = D * time;
% fprintf('N = %d\n', N);
% fprintf('Origin Integral %f\n', I);
% fprintf('Estimate of Integral %f\n', M);
% fprintf('Elapced time %f seconds\n', time);
% fprintf('Dispersion(*10^8) %f\n', D*10^8);
% fprintf('Laboriousness(*10^8) %f\n', S*10^8);
