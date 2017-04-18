clc
clear all
close all
N = 750; %simulation length
M = 5; %channel length
T = 20; %number of independent trials
cir = zeros(1,2*M-1); %cascade impulse response
MSE_vec = zeros(T,N-3); %Mean Square Error vector
MSE_f=zeros(1,N-3); %final MSE vector
o_p=zeros(1,N); %output from the equalizer in Tth trial
for j = 1 : T
u = sign(randn(1,N)); %input signal
c = randn(M,1); %channel to be equalized
c = c / norm(c);
z = filter(c,1,u); %channel output
SNR = 30; %Signal to Noise Ratio
var_v = var(z) * 10^(-SNR/10); %additive noise to the channel output
v = var_v^0.5 * randn(1,N);
y = z + v; %input to the equalizer
%-----------------------------LMS channel equalization-----------------------%
w = zeros(M,1);
y_d = zeros(1,M);
step =.12;
MSE=0;
for k = 4 : N
y_d = [y(k) y_d(1:M-1)];
e = u(k-3) - y_d * w;
w = w + step * y_d' * e ;
if j==T
o_p(1,k)=y_d * w;
end
MSE=MSE+(abs(e).^2);
MSE_avg=MSE/(k-3);
MSE_vec(j,:)=[MSE_avg MSE_vec(j,1:N-4)];%mse is calculated row wise,
                                        %MSE_vec(j,1:N-4)(j=row no,elements
                                        %are from 1st to N-4
end
MSE_vec(j,:)=fliplr(MSE_vec(j,:));%rows are inverted from left to right
q=[4:1:N];
figure(2);
plot(q,MSE_vec(j,:));
title('mean square error for different trial');
xlabel('No of iterations ');
ylabel('MSE');
hold on;
cir = cir + conv(w,c)';
end
hold off
figure(1);
subplot(2,1,1);
plot(q,MSE_vec(j,:)); 
title('MSE in Tth trial');
xlabel('No of iterations');
ylabel('MSE');
hold on;
o_p_b=sign(o_p); %converting to binary(+1/-1)
n=[1:N];
plot(n,v); %noise in the final trial
hold off;
temp=ones(T,1);
for j=1:N-3
MSE_f(1,j)=MSE_f(1,j)+((MSE_vec(:,j)'*temp)/T);
end
subplot(2,1,2);
plot(q,MSE_f);
title('MSE in final trial'); 
xlabel('No of iterations');
ylabel('MSE');
%-----------------------FOR SIGNALS IN FINAL TRIAL-------------------------%
figure(3);
subplot(3,1,1);
plot(c);
title('channel impulse response');
xlabel('time');
ylabel('magnitude');
subplot(3,1,2);
plot(w);
title('equalizer impulse response');
xlabel('time');
ylabel('magnitude');
subplot(3,1,3);
stem(cir/T);
title('cascade channel-equalizer impulse response');
ylabel('magnitude');
xlabel('taps');
figure(4);
subplot(2,2,1);
n_i=n(4:N);
stem(n_i,u(1:N-3));
xlabel('time');
title('input signal');
ylabel('magnitude');
subplot(2,2,2);
stem(n,z);
xlabel('time');
title('output after passing through channel');
ylabel('magnitude');
subplot(2,2,3);
n_o=n(1:N-3);
stem(n_o,o_p(1:N-3));
xlabel('time');
title('Final output');
ylabel('magnitude');
subplot(2,2,4);
stem(n_o,o_p_b(1:N-3));
xlabel('time');
title('Binary final output');
ylabel('magnitude');

