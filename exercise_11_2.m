%Exercise2
%Comparing the timing

clear
close all;


load('sampleEEGdata.mat');%Loading the EEG sampled data
eegdat4convol = EEG.data(47,:,1);% get data of one epoch from channel 47
 
% creating a Gaussian kernel
time = -1:1/EEG.srate:1;
s = 5/(2*pi*30);
gaussian = exp((-time.^2)/(2*s^2))/30;

tic
for ind = 1 : 1000
    
    %zero padding with length of gaussian kernel
    eegdat4conv_padded = [zeros(1,length(gaussian)-1) eegdat4convol zeros(1,length(gaussian)-1)];

    % used for cutting the result of convolution
    half_of_kernel_size = ceil((length(gaussian)-1)/2);
    % creating the output vector
    man_conv_result = zeros(1,length(eegdat4convol)+length(gaussian)-1);
    %flipping the gaussian kernel
    flip_g = fliplr(gaussian);

    for i = 1 : length(man_conv_result)
        man_conv_result(i) = sum(eegdat4conv_padded(i:i+length(gaussian)-1).*flip_g);
    end
    % cut off edges
    man_conv_result = man_conv_result(half_of_kernel_size+1:end-half_of_kernel_size);
end

telapsed_manual = toc;

save ('time_manual.mat', 'telapsed_manual');
%%
clear
close all;


load('sampleEEGdata.mat');%Loading the EEG sampled data
eegdat4convol = EEG.data(47,:,1);% get data of one epoch from channel 47
 
% creating a Gaussian kernel
time = -1:1/EEG.srate:1;
s = 5/(2*pi*30);
gaussian = exp((-time.^2)/(2*s^2))/30;

tic
for ind = 1 : 1000
    N_conv = length(eegdat4convol) + length(gaussian) - 1;

    fourier_eegdat4convol = zeros(1,N_conv);
    fourier_gaussian = zeros(1,N_conv);

    eegdat4convol_padded = [eegdat4convol, zeros(1, length(gaussian) -1)];
    gaussian_padded = [gaussian, zeros(1,length(eegdat4convol) -1)];
    half_of_kernel_size = ceil((length(gaussian)-1)/2);

    % frequency domain representation
    k = ((1:N_conv)-1)/N_conv;

    for fi = 1:N_conv
    
        sine_wave = exp(-1i*2*pi*(fi-1).* k);
        fourier_eegdat4convol (fi) = sum (sine_wave.*eegdat4convol_padded);
        fourier_gaussian (fi) = sum (sine_wave.*gaussian_padded);
    
    end

    fourier_conv = (fourier_eegdat4convol.* fourier_gaussian)/N_conv;

    % time domain representation
    time_conv_result = zeros(1,N_conv);

    for fi = 1:N_conv
        % scale sine wave by fourier coefficient
        sine_wave_conv = fourier_conv(fi)*exp(1i*2*pi*(fi-1).*k);
        % sum sine waves together (take only real part)
        time_conv_result = time_conv_result + real(sine_wave_conv);
    end
end

telapsed_DTFT = toc;

save ('time_DTFT.mat','telapsed_DTFT');
%% 
clear
close all;

load('sampleEEGdata.mat');%Loading the EEG sampled data
eegdat4convol = EEG.data(47,:,1);% get data of one epoch from channel 47
 
% creating a Gaussian kernel
time = -1:1/EEG.srate:1;
s = 5/(2*pi*30);
gaussian = exp((-time.^2)/(2*s^2))/30;


tic
for ind = 1:1000
    f_eegdat4convol = fft(eegdat4convol,length(eegdat4convol)+length(gaussian)-1); %FFT the data
    f_kernel = fft(gaussian,length(eegdat4convol)+length(gaussian)-1); %FFT the kernel
    fft_conv_result = real(ifft(f_eegdat4convol.*f_kernel)); %Multiply them together and perform IFFT
end

telapsed_FFT = toc;
save ('time_FFT.mat', 'telapsed_FFT');
%%
load time_manual.mat, 
load time_DTFT.mat
load time_FFT.mat

elapsed_time = [telapsed_manual, telapsed_DTFT, telapsed_FFT];
figure

bar(elapsed_time)