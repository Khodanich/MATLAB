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

Asin = 100;
Asin1 = 100;

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
phSignal = 0;
phSin = 10;
phSin1 = 50;

s1n_noise =  s1*exp(1i*phSignal/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');%% сигнал1 с шумом
s2n_noise =  s1*exp(1i*phSignal/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');%%сигнал 2 с шумом
sin_noise_1 =  SignalNoise*exp(1i*phSin/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');%% сигнал с шумом
sin_noise_2 =  SignalNoise1*exp(1i*phSin1/57.3) + wgn(1,lengthCorr,NoiseDb,'complex');%% сигнал с шумом

for j1_ctr = 1:1:360 
    
    Signal_sdvig_phase = s1n_noise*exp(1i*((deg2rad(j1_ctr - 1))));
    Signal_sdvig_phase_2 = Signal_sdvig_phase + sin_noise_1 + sin_noise_2;
    Energy(j1_ctr) = sum(Signal_sdvig_phase_2.^2);
    %Signal_sdvig_phase_21_sin = abs(Correlation(Signal_sdvig_phase_2,s)); 
    
    %q_phase_Quality_sin(j1_ctr) = Quality(Signal_sdvig_phase_21_sin,lengthCorr);
   
end

mass_temp = 0:359; 
plot(mass_temp,Energy);

