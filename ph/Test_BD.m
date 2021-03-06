close all;
system = 'gps';

sample = 1;
lengthCorr = 4092;

s1 = zeros(1,lengthCorr);
svnum1 = 2;
codeX = Code_Generator( system, sample, svnum1);%%������������ C/A ���
s = zeros(1,lengthCorr/2);%%�������� ���
r = 1;
for i = 1:1023 
   s(r) = codeX(i);
   r = r + 1;     
   s(r) = codeX(i);
   r = r + 1;
end

s1(1:lengthCorr/2) = s;
s1((lengthCorr/2 + 1):lengthCorr) = s;

NoiseDb = 25;
AmpSignal = 5;%%��������� �������
code1xTemp1 = s1;
code1xTemp = s1;

for i = 1:lengthCorr
    if s1(i) == 1
        s1(i) = AmpSignal;
        s(i) = AmpSignal;        
    else
        s1(i) = -AmpSignal;
        s(i) = -AmpSignal;
    end
end

Asin = 1000;
Asin1 = 1000;

Fd = 2.046*1000000;
Fd1 = 2.046*1000000;
Fs = 600*1000;
Fs1 = 600*1000;

SignalNoise = zeros(1,lengthCorr);
SignalNoise1 = zeros(1,lengthCorr);

for i = 1:lengthCorr        %%��������� �����
  SignalNoise(i) = Asin*sin(2*pi*Fs/Fd*i); 
  SignalNoise1(i) = Asin1*sin(2*pi*Fs1/Fd1*i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phref = 0;


phSignal1 = 50;
phSignal2 = 60;

phSin = 120;
phSin1 = 90;

s1_reference = s1*exp(1i*phref/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');                    %%������� ������
s1n_noise =  s1*exp(1i*phSignal1/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');%% ������1 � ����� ������� �� 30 ��������
%s2n_noise =  s1*exp(1i*phSignal2/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');%%������ 2 � ����� ������� �� 60 ��������



sin_noise_1 =  SignalNoise*exp(1i*phSin/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');%% ������1
%sin_noise_2 =  SignalNoise1*exp(1i*phSin1/57.3);% + wgn(1,lengthCorr,NoiseDb,'complex');%% ������2

%%% ������ ����� ���
s1_WhiteNoise = s1*exp(1i*phSignal1/57.3) + wgn(1,lengthCorr,NoiseDb,'complex'); 
%%���� ������
s_one_interference = s1*exp(1i*phSignal1/57.3) + SignalNoise*exp(1i*phSin/57.3) + wgn(1,lengthCorr,NoiseDb,'complex'); 
s1_reference_one_interf = s1*exp(1i*phref/57.3) +  SignalNoise + wgn(1,lengthCorr,NoiseDb,'complex');

%%��� ������

s1n_mix = s1*exp(1i*phSignal1/57.3) + SignalNoise*exp(1i*phSin/57.3)+ SignalNoise1*exp(1i*phSin1/57.3)+ wgn(1,lengthCorr,NoiseDb,'complex');     %%����� �������1 � �����
s2n_mix = s1*exp(1i*phSignal2/57.3) + SignalNoise*exp(1i*phSin/57.3)+ SignalNoise1*exp(1i*phSin1/57.3)+ wgn(1,lengthCorr,NoiseDb,'complex');   %%����� �������2 � �����
s_referense_interf = s1*exp(1i*phref/57.3) + SignalNoise + SignalNoise1 + wgn(1,lengthCorr,NoiseDb,'complex');%% ����� �������� ������� � �����


mass_temp = 0:359;
%%����������� ������ ������ ����
% for j1_ctr = 1:1:360     
%     s2n_noise_interf_WhiteNoise = s1_reference*exp(1i*((deg2rad(j1_ctr - 1))));
%     %Signal_sdvig_phase2 = s2n_noise*exp(1i*((deg2rad(j1_ctr - 1))));
%     s2n_noise_interf_sum = s2n_noise_interf_WhiteNoise + s1_WhiteNoise; 
%     %Signal_sdvig_phase_Sum = Signal_sdvig_phase2 + sin_noise_1 + sin_noise_2;
%     
%     Signal_sdvig_phase_22_sin = abs(Correlation(s2n_noise_interf_sum,s)); 
%   
%     %Energy1(j1_ctr) = sum(abs(s2n_noise_interf).^2);
%     q_phase_Quality_sin(j1_ctr) = Quality(Signal_sdvig_phase_22_sin,lengthCorr);   
%     %max_phase_max(j1_ctr) =  max(Signal_sdvig_phase_22_sin);
% end
 
%%����������� ����� ������
% for j1_ctr = 1:1:360     
%     s2n_noise_interf = s1_reference_one_interf*exp(1i*((deg2rad(j1_ctr - 1))));
%     %Signal_sdvig_phase2 = s2n_noise*exp(1i*((deg2rad(j1_ctr - 1))));
%     s2n_noise_interf_sum = s2n_noise_interf + s_one_interference; 
%     %Signal_sdvig_phase_Sum = Signal_sdvig_phase2 + sin_noise_1 + sin_noise_2;
%     
%     Signal_sdvig_phase_22_sin = abs(Correlation(s2n_noise_interf_sum,s)); 
%   
%     %Energy1(j1_ctr) = sum(abs(s2n_noise_interf).^2);
%     q_phase_Quality_sin(j1_ctr) = Quality(Signal_sdvig_phase_22_sin,lengthCorr);   
%     %max_phase_max(j1_ctr) =  max(Signal_sdvig_phase_22_sin);
% end
% 
% plot(mass_temp,q_phase_Quality_sin);

%%����������� ���� �����
for i = 1:1:360
    s1n_noise_interf = s1n_mix*exp(1i*((deg2rad(j - 1))));    
    for j = 1:1:360
        s2n_noise_interf = s2n_mix*exp(1i*((deg2rad(i - 1))));
        %Signal_sdvig_phase2 = s2n_noise*exp(1i*((deg2rad(j1_ctr - 1))));
        s2n_noise_interf = s1n_noise_interf + s2n_noise_interf + s_referense_interf; 
        %Signal_sdvig_phase_Sum = Signal_sdvig_phase2 + sin_noise_1 + sin_noise_2;

        %Signal_sdvig_phase_22_sin = abs(Correlation(s2n_noise_interf,s)); 
        Energy1(i) = sum(abs(s2n_noise_interf).^2);
        %q_phase_Quality_sin(j1_ctr) = Quality(Signal_sdvig_phase_21_sin,lengthCorr);    
    end  
    MaxEnergy(j) = max(Energy1);
end
[X,Y]=meshgrid(mass_temp,mass_temp);
plot3(X,Y,MaxEnergy);
% 
% plot(mass_temp,Energy1);
% figure;
% 
% for j1_ctr = 1:1:360 
%     
%     s1n_noise_interf = s1n_mix*exp(1i*((deg2rad(j1_ctr - 1))));
%     %Signal_sdvig_phase2 = s2n_noise*exp(1i*((deg2rad(j1_ctr - 1))));
%     %s1n_noise_interf = s1n_noise_interf + s_reference_interf; 
%     Signal_sdvig_phase_Sum = s1n_noise_interf + s_reference_interf;
%     
%     Signal_sdvig_phase_21_sin = abs(Correlation(Signal_sdvig_phase_Sum,s)); 
%     
%     %Energy(j1_ctr) = sum(abs(s1n_noise_interf).^2);
%     q_phase_Quality_sin(j1_ctr) = Quality(Signal_sdvig_phase_21_sin,lengthCorr);
%    
% end


% plot(mass_temp,q_phase_Quality_sin);

