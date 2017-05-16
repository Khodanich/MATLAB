close all;
system = 'gps';

sample = 1;
lengthCorr = 4092;

s1 = zeros(1,lengthCorr);
svnum1 = 2;
codeX = Code_Generator( system, sample, svnum1);%%������������ C/A ���
s = zeros(1,lengthCorr/2);%%�������� ���

m = mod(13,360);
deg2rad(m);
% svnum2 = 1;
% code2 = Code_Generator( system, sample, svnum2);
r = 1;
for i = 1:1023 
   s(r) = codeX(i);
   r = r + 1;     
   s(r) = codeX(i);
   r = r + 1;
end

plot(Correlation(s, s));
ylabel('���������');
xlabel('�����');
grid on;
%title('������������������ ������� ��������� �������');% ������� �������
s1(1:lengthCorr/2) = s;
s1((lengthCorr/2 + 1):lengthCorr) = s;%%��������  ������ � �������� �������������� � ��� ���� ������
% code1x((lengthCorr/2 + 1):lengthCorr*3/4) = code1;
% code1x((lengthCorr/4*3 +1):lengthCorr) = code1;
Asin = 1000;%% ��������� ������
AmpSignal = 1;%%��������� �������
NoiseDb = 25; %%�������� ������

for i = 1:lengthCorr
    if s1(i) == 1
        s1(i) = AmpSignal;
        s(i) = AmpSignal;        
    else
        s1(i) = -AmpSignal;
        s(i) = -AmpSignal;
    end
end

%%% s �������� ������ 
code1xTemp1 = s1;
code1xTemp = s1;


z = 1:lengthCorr;


cor11 = abs(Correlation(s1, s));
max(cor11);
qq = fft(s);
grapSilent = cor11;
%plot(fftshift(abs(cor11)));

plot (grapSilent);
axis([0,4080,0,4500]);
ylabel('���������');
xlabel('�������');
%title('�������������� ������� �������� s � s1');% ������� �������
grid on;
figure;

%% ������ ������

Fd = 2.046*1000000;
Fs = 500*1000;
Fs1 = 600*1000;
FftL = 4092;

SignalNoise = zeros(1,lengthCorr);
for i = 1:lengthCorr
  SignalNoise(i) = Asin*sin(2*pi*Fs/Fd*i); %%�������� ������ 
  if Asin == 0
   SignalNoise(i) = 0; 
  end
end

%%%
% for i = 1:lengthCorr
%   SignalNoise15(i) = Asin*sin(2*pi*Fs1/Fd*i); %%�������� ������ 
%   if Asin == 0
%    SignalNoise15(i) = 0; 
%   end
% end

%%%
SperctrSin = abs(fft(SignalNoise,FftL));

F = 0:Fd/FftL:Fd/2-1/FftL;
plot(F,SperctrSin(1:length(F)));
title('������ ������');
axis([0,1000000,0,2500000]);
grid on;
xlabel('�������, ��');
ylabel('���������');

%% ��������� ������
% code1xTemp = abs(fft(code1xTemp));
for i = 1:lengthCorr
    code1xTemp(i) = SignalNoise(i) + code1xTemp(i);% + SignalNoise15(i);
end
figure;
temp11 = abs(fft(code1xTemp,FftL)); 
plot(F,temp11(1:length(F)));
% plot(abs(fft(code1xTemp)));
% plot(F,code1xTemp(1:length(F)));
xlabel('�������, ��');
ylabel('���������');
grid on;
title('������ ������-������(������� ��������� 500 ���)');% ������� �������
figure;
%%
corrr = Correlation(code1xTemp,s);
corrr = abs(corrr);
plot(corrr);
axis([0,4080,0,70000]);
ylabel('���������');
xlabel('�������');
%title('���������� ��������� ������� � �������� � ������� ');% ������� �������
grid on;
%% ���������� �������� q ����� ���� ���
qNoise = Quality(corrr,lengthCorr);
%%
%%���������� �������� q ��� ����
corSilent  = abs(cor11);

figure;
plot (corSilent);
title('���������� ��������� ������� ��� ������ ');% ������� �������
qSilent = Quality(corSilent,lengthCorr);
%% ����������� �������� �� ������ �������
 SignalNoise_temp = SignalNoise/Asin;
AsinAm = zeros(1,lengthCorr);

for j = 1:lengthCorr    
    SignalNoise1 = SignalNoise_temp*(j*2);
   % code1xTemp1 = (wgn(1,lengthCorr,ctr1/10)) + (code1xTemp1);
    code1xTemp1 = (SignalNoise1) + (code1xTemp1);
    SignalNoise1 = abs(Correlation(code1xTemp1,s));
    AsinAm(j) = (j*2);%%��������� ������    
    qNoise1(j) = Quality(SignalNoise1,lengthCorr);
end

ASingnalAm = ones(1,lengthCorr);
ASingnalAm = ASingnalAm*AmpSignal;
for  i3 = 1:lengthCorr
    SignalNoise1(i3) = 20*log(ASingnalAm(i3)/AsinAm(i3));
end;
%%
figure;
% plot (qNoise1,SignalNoise1);
plot (SignalNoise1,qNoise1);
ylabel('���������� ��������');
xlabel('��������� ������/���, ��');
%title('������ ����������� �������� ������� �� ��������� ������/��� ');% ������� �������
grid on;
figure;
%% ������� ����
j11 = 50;
j22 = 0;
j11_sin = 120;
% for j1 = 1:lengthCorr
%     z_code1x(j1) = z_code1x(j1)*exp(1i*j11);
%     z_code1x_real(j1) = real(z_code1x(j1));
%     z_code1x_imag(j1) = imag(z_code1x(j1));
%     phase_main(j1) = atan2(z_code1x_imag(j1),z_code1x_real(j1))*57.3; 
% end

% y11 = code1x ;
% y22 = (code1x)*exp(1i*180/57.3);
% z33 = y11 + y22;
% z22 = abs(Correlation(z33,y11));
% plot(z22);
% q_ttr = Quality(z22,lengthCorr);
%% �������� �������� � ����
mass_temp = 0:359; 

z231 = wgn(1,lengthCorr,25);
z233 = awgn(s1,25,'measured');
%Signal_sdvig_1 =  code1x + wgn(1,lengthCorr,10);   %������ ������ ��������� � ���� 0, ��� ������ ������� ����� ������� wgn(1,lengthCorr,25);

%s1n =  s1*exp(1i*j22/57.3) + wgn(1,lengthCorr,NoiseDb);
%s1n_cl = s1*exp(1i*j11_sin/57.3);
s1n =  s1*exp(1i*j22/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');
re1 = abs(Correlation(s1n,s));
plot(re1);
axis([0,4080,0,4000]);
grid on;
ylabel('���������');
xlabel('�������');
%title('�������������� ������� ��������� ������� � ������ ������� � ������ ����');% ������� �������

re2 = Quality(re1,lengthCorr);
% q_s1n = Quality(abs(Correlation(s1n,s)),lengthCorr);
% q_s1n_sin = Quality(abs(Correlation(s1n_sin,s)),lengthCorr);
%s1n =  s1*exp(1i*j22/57.3) + SignalNoise;
%s1n = awgn(s1,25,'measured');
%%%
%s2np = s1*exp(1i*j11/57.3) + wgn(1,lengthCorr,NoiseDb);%%������ ������ ��������� � ����� �� ����, ������ ������� ��� � ����
%s2np_cl = s1*exp(1i*j11/57.3);
%s2np = s1*exp(1i*j11/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');
s2np_sin = s1*exp(1i*j11/57.3) + wgn(1,lengthCorr,NoiseDb,'complex') + SignalNoise*exp(1i*j11_sin/57.3);
s1n_sin = s1*exp(1i*j22/57.3) + wgn(1,lengthCorr,NoiseDb,'complex') + SignalNoise;
% q_s2np = Quality(abs(Correlation(s2np,s)),lengthCorr);
% q_s2np_sin = Quality(abs(Correlation(s2np_sin,s)),lengthCorr);
%s2np = s1*exp(1i*j11/57.3) + SignalNoise;
%s2np = awgn(s1,25,'measured');
% q_phase_Quality = zeros(1,360);
% 
figure;
plot(abs(Correlation(s1n_sin, s))); 
% figure;
% s1n1 =  s1n*exp(1i*30/57.3) + s2np;
% re1 = abs(Correlation(s1n1,s)); 
% plot(re1);
% qww = Quality(re1,lengthCorr);
% title('�������� �������� � ����');% ������� �������
% 
% figure;
% s2n2 =  s1n*exp(1i*210/57.3)  + s2np;
% re2 = abs(Correlation(s2n2,s)); 
% plot(re2);
% qww1 = Quality(re2,lengthCorr);
% title('�������� �������� � ����������� ');% ������� �������
% figure;
% t = 1;
% S_temp = s1n + s2np;
% re = abs(Correlation(S_temp,s1));
% qww = Quality(re,lengthCorr);
% plot(re);
% title('���������� ����� �������� s1n � s2np � �������� �������� s1');% ������� �������

% figure;
% plot(abs(Correlation(Signal_sdvig_1,code1x)));
% title('���������� ��������� ������� � �������� s1n ');% ������� �������
% 
% qw2 = Quality(abs(Correlation(Signal_sdvig_1,code1x)),lengthCorr);
% 
% figure;
% plot(abs(Correlation(Signal_stay_2_phase,code1x)));
% title('���������� ��������� ������� � �������� s2np');% ������� �������
% 
% qw21 = Quality(abs(Correlation(Signal_stay_2_phase,code1x)),lengthCorr);

% % % ttt = abs(Correlation(Signal_sdvig_1,code1x));
% % % plot(ttt);
% % % qou_1 = Quality(ttt,lengthCorr);
% % % title('���������� ��������� ������� s � s1p  ');
% % % figure;
% % % ttt = abs(Correlation(Signal_stay_2_phase,code1x));
% % % plot(ttt);
% % % qou_2 = Quality(ttt,lengthCorr);
% % % 
% % % title('���������� ��������� ������� s � s2p ���������� �� ���� �� 30 ��������  ');
Signal_sdvig_phase = zeros(1,lengthCorr);
Signal_sdvig_phase_2 = zeros(1,lengthCorr);
Signal_sdvig_phase_21_sin = zeros(1,lengthCorr);
% for i22 = 1:4
%     s1n =  s1*exp(1i*j22/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');
for j1_ctr = 1:1:360    
  %%  wgn + sin 
  Signal_sdvig_phase = s1n_sin*exp(1i*((deg2rad(j1_ctr - 1))));
  Signal_sdvig_phase_2 = Signal_sdvig_phase + s2np_sin;
  Signal_sdvig_phase_21_sin = abs(Correlation(Signal_sdvig_phase_2,s));  
  q_phase_Quality_sin(j1_ctr) = Quality(Signal_sdvig_phase_21_sin,lengthCorr);
  %% only wgn
%   Signal_sdvig_phase = s1n*exp(1i*((deg2rad(j1_ctr - 1))));
%   Signal_sdvig_phase_2 = Signal_sdvig_phase + s2np;
%   Signal_sdvig_phase_21 = abs(Correlation(Signal_sdvig_phase_2,s));  
%   q_phase_Quality(j1_ctr) = Quality(Signal_sdvig_phase_21,lengthCorr);  
%   %% only signal
%    Signal_sdvig_phase = s1n_cl*exp(1i*((deg2rad(j1_ctr - 1))));
%   Signal_sdvig_phase_2_cl = Signal_sdvig_phase + s2np_cl;
%   Signal_sdvig_phase_21_cl = abs(Correlation(Signal_sdvig_phase_2_cl,s));  
%   q_phase_Quality_cl(j1_ctr) = Quality(Signal_sdvig_phase_21_cl,lengthCorr);
%   
  
  %%
  

%  
%   if j1_ctr == 210
%      t = 1; 
%      d3 = mean(Signal_sdvig_phase_2);
%      plot(abs(Correlation(Signal_sdvig_phase_2,s1)));
%        axis([0,4100,0,4900]);
%      ylabel('���������');
%      xlabel('�������');
%      grid on;
%     % title('���������� ��������� ������� � ������� ���������� �� ���� �� 210 ��������');% ������� �������
%   end
% if j1_ctr == 10  
%     
%     
% end
% 
% if j1_ctr == 30    
%     figure;
%     z123 = abs(Correlation(Signal_sdvig_phase_2_cl,s));  
%     q_phase_Quality_cl(j1_ctr) = Quality(Signal_sdvig_phase_21_cl,lengthCorr);
%     plot(z123);
%     grid on;
%     title('���������� ��������� ������� � ������� ���������� �� ���� �� 30 ��������');% ������� �������
% 
% end
% % % 
% if j1_ctr == 210 
%     figure;
%     plot(abs(Correlation(Signal_sdvig_phase_2,s)));  
%     axis([0,4100,0,4900]);
% %    q_phase_Quality(j1_ctr) = Quality(Signal_sdvig_phase_21,lengthCorr);
%     grid on;
%     
%     %title('�������������� ������� �������� s1 � s2np � �����������');% ������� �������
% end

% if j1_ctr == 50    
%     figure;
%     plot(Signal_sdvig_corr);
%     title('150');% ������� �������
% end
% 
% 
% if j1_ctr == 210  
%     figure;
%     plot(Signal_sdvig_corr);
%     t1 = max(Signal_sdvig_corr);
%     zzz1 = mean(Signal_sdvig_corr);
%     z = mean(Signal_sdvig_corr) - (max(Signal_sdvig_corr) + Signal_sdvig_corr(2))/lengthCorr;
%     title('330 ');% ������� �������
% end


%%


 
 max_phase_max(j1_ctr) =  max(Signal_sdvig_phase_21_sin);
%  max_phase = max(Signal_sdvig_phase_21_sin);
%  
%   for i = 1:lengthCorr
%     if  Signal_sdvig_phase_21_sin(i) == max_phase
%         max_Num(j1_ctr) = i;
%     end
%   end
end
%%
% % % 
% % % for ctr1 = 1:500
% % %    signal11 = code1x + wgn(1,lengthCorr,ctr1/10); 
% % %    a_gauss(ctr1) = ctr1/10;
% % %    Z1 = abs(Correlation(signal11,(code1x)));
% % %    q_temp_22(ctr1) = Quality(Z1,lengthCorr);
% % % end
% % % end
% % % 
% % %     for  i3 = 1:500
% % %         SignalNoise_1(i3) = 20*log(AmpSignal)- a_gauss(i3);  
% % %     end;  
% % %     
% % % plot(SignalNoise_1,q_temp_22);
%%
% plot(mass_temp,q_phase_Quality_cl);
% ylabel('�������� ����������');
% xlabel('����');
% grid on;
% title('����������� ��������� ���������� �� ���� ');% ������� �������
% figure;
%%
% plot(mass_temp,max_phase_max);
% ylabel('�������� ����������');
% xlabel('����');
% grid on;
% title('����������� ��������� ���������� �� ���� ');% ������� �������
% figure;
% figure;
% plot(mass_temp,max_phase_max);
% ylabel('��������');
% xlabel('����');
% title('����������� ��������� ���������� ������� �� ���� ');% ������� �������
% %%%
% subplot(2, 2, i22),plot(mass_temp,q_phase_Quality);
% 
% plot(mass_temp,q_phase_Quality);
% axis([0,370,60,100]);
% ylabel('��������');
% xlabel('����');
% grid on;
% %title('����������� �������� ������� �� ���� ��� ����������� ������ ���� ');% ������� �������
% end;
% figure;
%%%%


polar(mass_temp*pi/180,q_phase_Quality_sin/100);

figure;
plot(mass_temp,q_phase_Quality_sin);
%axis([0,370,35,100]);
% ylabel('��������');   
ylabel('��������');
xlabel('����');
grid on;
title('����������� �������� ������� �� ���� ��� ����������� ������ ���� � ������');% ������� �������
% code1x_temp = code1x + code1x;
% code1x_phase = code1x*exp(1i*180/57.3);
% code1x_temp = code1x + code1x_phase;
% phase = atan2(imag(code1x_temp),real(code1x_temp))*57.3;

%% 

    

