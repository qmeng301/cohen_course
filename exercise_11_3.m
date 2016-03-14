%Exercise3
clear
close all
clc

f1 = 10;
f2 = 20;
f3 = 30;

fs = 1000;
T = 10;

t = (1:fs*T)/fs; %defining the time points
f_res = fs/length(t);
f_axis = (1:length(t))*f_res;
s1 = sin(2*pi*f1*t);
s2 = sin(2*pi*f2*t);
s3 = sin(2*pi*f3*t);

sine_waves = s1+s2+s3; 

figure(1)
subplot(3,1,1)
plot(t,s1);
subplot(3,1,2)
plot(t,s2);
subplot(3,1,3)
plot(t,s3);

%%
fft_sine_waves = abs(fft(sine_waves));
figure(2)
subplot(2,1,1)
plot(t, sine_waves);
subplot(2,1,2)
plot(f_axis(1:length(fft_sine_waves)/2+1),fft_sine_waves(1:length(fft_sine_waves)/2 + 1));
%%
noise_small = 0.5*randn(1,length(t));
sine_waves_noise_s = noise_small + sine_waves;
fft_sine_waves_noise_s = abs(fft(sine_waves_noise_s));

figure(3)
subplot(2,1,1)
plot(t, sine_waves_noise_s);
subplot(2,1,2)
plot(f_axis(1:length(fft_sine_waves)/2+1),fft_sine_waves_noise_s(1:length(fft_sine_waves)/2 + 1));
%%
noise_large = 5*randn(1,length(t));
sine_waves_noise_l = noise_large + sine_waves;
fft_sine_waves_noise_l = abs(fft(sine_waves_noise_l));

figure(4)
subplot(2,1,1)
plot(t, sine_waves_noise_l);
subplot(2,1,2)
plot(f_axis(1:length(fft_sine_waves)/2+1),fft_sine_waves_noise_l(1:length(fft_sine_waves)/2 + 1));