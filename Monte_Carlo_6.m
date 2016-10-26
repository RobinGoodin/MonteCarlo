clc
clear
close
%% Задание основных параметров
set(0,'RecursionLimit',10000)
rand('seed',666);
N = 300;
N_start  = 10;
N_step = 10;
C = 16.1975;
Ni = N_start:N_step:N;
res_old = zeros(length(Ni),1);
res_new = zeros(length(Ni),1);

func = @(x,y,z,u) exp(x*y*z*u)/sqrt(x*y*z*u);
f = @(x,y,z,u) 16*exp(x*y*z*u);
p = @(x,y,z,u) 1/(16*sqrt(x*y*z*u));
I = 0;
for k=0:1000
   I = I + 16*1/(factorial(k)*(2*k+1)^4);
end
y = rand(N,4);
%% Вычисление последовательности Холтона
p2 = zeros(1, N);
p3 = zeros(1, N);
p5 = zeros(1, N);
for k = 1:N
    M = k + 2;
    a2 = zeros(1, M);
    a3 = zeros(1, M);
    a5 = zeros(1, M);
    
    a2(1) = mod(k, 2);
    p2(k) = a2(1)/2;
    
    a3(1) = mod(k, 3);
    p3(k) = a3(1)/3;

    a5(1) = mod(k, 5);
    p5(k) = a5(1)/5;
    
    for j = 2:M
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

%%
for k = 1:length(Ni)
    
    % Старый алгоритм
    M_old = 0;
    for i=1:Ni(k)
        M_old = M_old + C + 16*(exp(y(i,1)^2*y(i,2)^2*y(i,3)^2*y(i,4)^2) ...
                - 1 - y(i,1)^2*y(i,2)^2*y(i,3)^2*y(i,4)^2) ;
    end
    M_old = M_old / Ni(k);
    res_old(k) = log(abs(M_old - I))/log(10);
    
    % Новый алгоритм
    M_new = 0;
    for i=1:Ni(k)
        M_new = M_new + C + 16*(exp((i/Ni(k))^2*p2(i)^2*p3(i)^2*p5(i)^2) ...
                - 1 - (i/Ni(k))^2*p2(i)^2*p3(i)^2*p5(i)^2) ;
    end
    M_new = M_new / Ni(k);
    res_new(k) = log(abs(M_new - I))/log(10);
    
end

%% Построени графиков
grid on;
hold on
    plot(log(Ni)/log(10), res_old, 'DisplayName', 'old');
    plot(log(Ni)/log(10), res_new, 'DisplayName', 'new');
hold off
xlabel('ln N');
ylabel('ln|teta - I|');
legend('show');
title('Plot');
