%
% Takes a list of hex numbers and displays an image
%
% Save the floating point STFT data as fftdata.dat in code composer studio,
% whilst paused in debug mode to view STFT.
%
% J.Harris, 14-3-2014
%
close all
clc
clear stft_image

N = 80000; % XLEN in C code, Length of the time domain signal
window_size = 1024;
freq_bins = 1024;
overlap = 0;
time_step = (window_size-overlap);
%time_blocks = (N-overlap) / (time_step);%freq_bins;
time_blocks = 72;%36;%floor((N-freq_bins)/(window_size-overlap)) + 1;


step = 2; % Step between beginning of the real part of the next number
a=1;


%%% for loop below for 50% overlap
for n=1:time_blocks
    offset = ((n-1)*2*freq_bins/2) + 1;% To calculate where the 'latest' fft is
    stft_image(1:512,n) = fftdata(offset:2:offset+1023) + (fftdata(offset+1:2:offset+1024) * j);
end

% for m=1:(length(fftdata)/2-1)
%     audio(m) = fftdata(2*m) + j * fftdata(2*m + 1);        
% end

% for n=1:time_blocks
%     offset = ((n-1)*2*freq_bins) + 1;% To calculate where the 'latest' fft is
%     stft_image(1:1024,n) = fftdata(offset:2:offset+2047) + (fftdata(offset+1:2:offset+2048) * j);
% end

stft_image = abs(stft_image);
imagesc(20*log(stft_image));
disp('Done!');