clear
clc
close
rand('seed',666);
%% ?????? ?????????
I = 0;
X = [0 0 0 0; 0.5 0.5 0.5 0.5; 1 1 1 1];
for k=0:100
   Xx = X(3,1)^(k+1/2)-X(1,1)^(k+1/2);
   Xy = X(3,2)^(k+1/2)-X(1,2)^(k+1/2);
   Xz = X(3,3)^(k+1/2)-X(1,3)^(k+1/2);
   Xu = X(3,4)^(k+1/2)-X(1,4)^(k+1/2);
   I = I + 16*(Xx*Xy*Xz*Xu)/(factorial(k)*(2*k+1)^4);
end
vpa(I,10);
%I = 16.211553

%% ???????????? ???????? ????????? ??????????
I01 = [0.0 0.0 0.0 0.0; 0.5 0.5 0.5 0.5];
I02 = [0.0 0.0 0.0 0.5; 0.5 0.5 0.5 1.0];
I03 = [0.0 0.0 0.5 0.0; 0.5 0.5 1.0 0.5];
I04 = [0.0 0.0 0.5 0.5; 0.5 0.5 1.0 1.0];
I05 = [0.0 0.5 0.0 0.0; 0.5 1.0 0.5 0.5];
I06 = [0.0 0.5 0.0 0.5; 0.5 1.0 0.5 1.0];
I07 = [0.0 0.5 0.5 0.0; 0.5 1.0 1.0 0.5];
I08 = [0.0 0.5 0.5 0.5; 0.5 1.0 1.0 1.0];
I09 = [0.5 0.0 0.0 0.0; 1.0 0.5 0.5 0.5];
I10 = [0.5 0.0 0.0 0.5; 1.0 0.5 0.5 1.0];
I11 = [0.5 0.0 0.5 0.0; 1.0 0.5 1.0 0.5];
I12 = [0.5 0.0 0.5 0.5; 1.0 0.5 1.0 1.0];
I13 = [0.5 0.5 0.0 0.0; 1.0 1.0 0.5 0.5];
I14 = [0.5 0.5 0.0 0.5; 1.0 1.0 0.5 1.0];
I15 = [0.5 0.5 0.5 0.0; 1.0 1.0 1.0 0.5];
I16 = [0.5 0.5 0.5 0.5; 1.0 1.0 1.0 1.0];
II = I01;
II(:,:,2) = I02;
II(:,:,3) = I03;
II(:,:,4) = I04;
II(:,:,5) = I05;
II(:,:,6) = I06;
II(:,:,7) = I07;
II(:,:,8) = I08;
II(:,:,9) = I09;
II(:,:,10) = I10;
II(:,:,11) = I11;
II(:,:,12) = I12;
II(:,:,13) = I13;
II(:,:,14) = I14;
II(:,:,15) = I15;
II(:,:,16) = I16;

%% ?????????? ????????? ?????????? ????????????
I2 = zeros(1,16);
for i=1:16
    for k=0:100
        Xx = II(2,1,i)^(k+1/2)-II(1,1,i)^(k+1/2);
        Xy = II(2,2,i)^(k+1/2)-II(1,2,i)^(k+1/2);
        Xz = II(2,3,i)^(k+1/2)-II(1,3,i)^(k+1/2);
        Xu = II(2,4,i)^(k+1/2)-II(1,4,i)^(k+1/2);
        I2(1,i) = I2(1,i) + 16*(Xx*Xy*Xz*Xu)/(factorial(k)*(2*k+1)^4);
    end
end
p = zeros(1,16);
for i=1:16
    Xx = II(2,1,i)^(1/2)-II(1,1,i)^(1/2);
    Xy = II(2,2,i)^(1/2)-II(1,2,i)^(1/2);
    Xz = II(2,3,i)^(1/2)-II(1,3,i)^(1/2);
    Xu = II(2,4,i)^(1/2)-II(1,4,i)^(1/2);
    p(i) = Xx*Xy*Xz*Xu;
end
summa = sum(p);

%% ?????????? ??????
tic
NN = 100000;
N = zeros(1,16);
for i=1:16
    N(i) = round(NN*p(i));
end
M = zeros(1,16);
M2 = zeros(1,16);
D = zeros(1,16);
for k=1:16
    g = rand(N(k),4);
    for i=1:N(k)
        x = (II(1,1,k)^(1/2) + g(i,1)*(II(2,1,k)^(1/2)-II(1,1,k)^(1/2)))^2;
        y = (II(1,2,k)^(1/2) + g(i,2)*(II(2,2,k)^(1/2)-II(1,2,k)^(1/2)))^2;
        z = (II(1,3,k)^(1/2) + g(i,3)*(II(2,3,k)^(1/2)-II(1,3,k)^(1/2)))^2;
        u = (II(1,4,k)^(1/2) + g(i,4)*(II(2,4,k)^(1/2)-II(1,4,k)^(1/2)))^2;
        M(k) = M(k) + p(k)*16*exp(x*y*z*u);
        M2(k) = M2(k) + (p(k)*16*exp(x*y*z*u))^2;
    end
    M(k) = M(k) / N(k);
    M2(k) = M2(k) / N(k);
    D(k) = M2(k) - M(k)^2;
end
M_I = sum(M);
D_I = 0;
for j = 1:16
    D_I = D_I + D(j)*p(j)*p(j)/N(j);
end
time = toc;
S_I = D_I * time;
fprintf('N = %d\n', NN);
fprintf('Origin Integral %f\n', I);
fprintf('Estimate of Integral %f\n', M_I);
fprintf('Elapced time %f seconds\n', time);
fprintf('Dispersion(*10^8) %f\n', D_I*10^8);
fprintf('Laboriousness(*10^8) %f\n', S_I*10^8);

tic
GG = 0;
for k = 1:16
    GG = GG + sqrt(D(k))'*p(k);
end
N = zeros(1,16);
for i=1:16
    N(i) = round(NN*p(i)*sqrt(D(i))/(GG));
end
M = zeros(1,16);
M2 = zeros(1,16);
D = zeros(1,16);
for k=1:16
    g = rand(N(k),4);
    for i=1:N(k)
        x = (II(1,1,k)^(1/2) + g(i,1)*(II(2,1,k)^(1/2)-II(1,1,k)^(1/2)))^2;
        y = (II(1,2,k)^(1/2) + g(i,2)*(II(2,2,k)^(1/2)-II(1,2,k)^(1/2)))^2;
        z = (II(1,3,k)^(1/2) + g(i,3)*(II(2,3,k)^(1/2)-II(1,3,k)^(1/2)))^2;
        u = (II(1,4,k)^(1/2) + g(i,4)*(II(2,4,k)^(1/2)-II(1,4,k)^(1/2)))^2;
        M(k) = M(k) + p(k)*16*exp(x*y*z*u);
        M2(k) = M2(k) + (p(k)*16*exp(x*y*z*u))^2;
    end
    M(k) = M(k) / N(k);
    M2(k) = M2(k) / N(k);
    D(k) = M2(k) - M(k)^2;
end
M_I = sum(M);
D_I = 0;
for j = 1:16
    D_I = D_I + D(j)*p(j)*p(j)/N(j);
end

time = toc;
S_I = D_I * time;
fprintf('N = %d\n', NN);
fprintf('Origin Integral %f\n', I);
fprintf('Estimate of Integral %f\n', M_I);
fprintf('Elapced time %f seconds\n', time);
fprintf('Dispersion(*10^8) %f\n', D_I*10^8);
fprintf('Laboriousness(*10^8) %f\n', S_I*10^8);