clear
close all
clc

%%
%10.1
% inverted u
kernel_1 = [0.12 0.4 0.66 0.84 1 0.84 0.66 0.4 0.12]; % 9 points

% decay function
kernel_2 = [1 0.78 0.65 0.44 0.32 0.2 0.15]; % 7 points

%%
%10.2
load('../sampleEEGdata.mat');%Loading the EEG sampled data
eegdat4convol = EEG.data(47,51:100,1);% get 50 time points from EEG data from channel 47

% zero padding with length of kernel_1
eegdat4conv_padded_1 = [zeros(1,length(kernel_1)-1) eegdat4convol zeros(1,length(kernel_1)-1) ];

% zero padding with length of kernel_2
eegdat4conv_padded_2 = [zeros(1,length(kernel_2)-1) eegdat4convol zeros(1,length(kernel_2)-1) ];

% used for cutting the result of convolution
half_of_kernel_size_1 = ceil((length(kernel_1)-1)/2); % 4 points

half_of_kernel_size_2 = ceil((length(kernel_2)-1)/2); % 3 points

% initialize convolution output
convolution_result_1 = zeros(1,length(eegdat4convol)+length(kernel_1)-1); % 58

convolution_result_2 = zeros(1,length(eegdat4convol)+length(kernel_2)-1); % 56

for i=1:length(convolution_result_1)
    convolution_result_1(i) = sum(eegdat4conv_padded_1(i:i+length(kernel_1)-1).*kernel_1(end:-1:1));
end


for j=1:length(convolution_result_2)
    convolution_result_2(j) = sum(eegdat4conv_padded_2(j:j+length(kernel_2)-1).*kernel_2(end:-1:1));
end

% cut off edges
convolution_result_1 = convolution_result_1(half_of_kernel_size_1+1:end-half_of_kernel_size_1); % 50 points
convolution_result_2 = convolution_result_2(half_of_kernel_size_2+1:end-half_of_kernel_size_2); % 50 points

figure(1)

subplot(3,1,1)
plot(kernel_1,'-o');
title ('Kernel 1: inverted U')

subplot(3,1,2)
plot(eegdat4convol);
title ('EEG.data(47,50:100,1)')

subplot(3,1,3)
plot(convolution_result_1./sum(kernel_1))
title ('convolution result')


figure(2)

subplot(3,1,1)
plot(kernel_2,'-o');
title ('Kernel 2: Decay function')

subplot(3,1,2)
plot(eegdat4convol);
title ('EEG.data(47,50:100,1)')

subplot(3,1,3)
plot(convolution_result_2./sum(kernel_2))
title ('convolution result') 