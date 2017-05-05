close all;
system = 'gps';

sample = 1;
lengthCorr = 4092;

s1 = zeros(1,lengthCorr);
svnum1 = 2;
codeX = Code_Generator( system, sample, svnum1);%%генерируемый C/A код
s = zeros(1,lengthCorr/2);%%исходный код
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
AmpSignal = 10;%%амплитуда сигнала
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

for i = 1:lengthCorr        %%генерация помех
  SignalNoise(i) = Asin*sin(2*pi*Fs/Fd*i); 
  SignalNoise1(i) = Asin1*sin(2*pi*Fs1/Fd1*i);
end

for i = 1:lengthCorr
    code1xTemp(i) = SignalNoise(i) +  SignalNoise1(i) + code1xTemp(i) ;
end

corrr = Correlation(code1xTemp,s);
corrr = abs(corrr);
plot(corrr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phref = 180;
phSignal1 = 30;
phSignal2 = 60;
phSin = 30;
phSin1 = 80;

s1_reference = s1 + wgn(1,lengthCorr,NoiseDb,'complex');                    %%опорный сигнал
s1n_noise =  s1*exp(1i*phSignal1/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');%% сигнал1 с шумом сдвинут на 30 градусов
s2n_noise =  s1*exp(1i*phSignal2/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');%%сигнал 2 с шумом сдвинут на 60 градусов

sin_noise_1 =  SignalNoise*exp(1i*phSin/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');%% помеха1
sin_noise_2 =  SignalNoise1*exp(1i*phSin1/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');%% помеха2

s1n_mix = s1n_noise + sin_noise_1 + sin_noise_2;    %%сумма сигнала1 и помех
s2n_mix = s2n_noise + sin_noise_1 + sin_noise_2;    %%сумма сигнала2 и помех
s_reference_interf = s1_reference + sin_noise_1 + sin_noise_2;%% сумма опорного сигнала и помех

mass_temp = 0:359;

for j1_ctr = 1:1:360     
    s2n_noise_interf = s2n_mix*exp(1i*((deg2rad(j1_ctr - 1))));
    %Signal_sdvig_phase2 = s2n_noise*exp(1i*((deg2rad(j1_ctr - 1))));
    s2n_noise_interf = s2n_noise_interf + s_reference_interf; 
    %Signal_sdvig_phase_Sum = Signal_sdvig_phase2 + sin_noise_1 + sin_noise_2;
    
    %Signal_sdvig_phase_22_sin = abs(Correlation(s2n_noise_interf,s)); 
    Energy1(j1_ctr) = sum(abs(s2n_noise_interf).^2);
    %q_phase_Quality_sin(j1_ctr) = Quality(Signal_sdvig_phase_21_sin,lengthCorr);   
end


plot(mass_temp,Energy1);
figure;

for j1_ctr = 1:1:360 
    
    s1n_noise_interf = s1n_mix*exp(1i*((deg2rad(j1_ctr - 1))));
    %Signal_sdvig_phase2 = s2n_noise*exp(1i*((deg2rad(j1_ctr - 1))));
    %s1n_noise_interf = s1n_noise_interf + s_reference_interf; 
    Signal_sdvig_phase_Sum = s1n_noise_interf + s_reference_interf;
    
    Signal_sdvig_phase_21_sin = abs(Correlation(Signal_sdvig_phase_Sum,s)); 
    
    %Energy(j1_ctr) = sum(abs(s1n_noise_interf).^2);
    q_phase_Quality_sin(j1_ctr) = Quality(Signal_sdvig_phase_21_sin,lengthCorr);
   
end


plot(mass_temp,q_phase_Quality_sin);

