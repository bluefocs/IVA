%
% Takes a list of hex numbers and displays an image
%
%
close all
clc
clear stft_image

N = 80000; % XLEN in C code, Length of the time domain signal
window_size = 1024;
freq_bins = (window_size / 2 ) + 1; % 
overlap = 0;
time_step = (window_size-overlap);
%time_blocks = (N-overlap) / (time_step);%freq_bins;
time_blocks = floor((N-freq_bins)/(window_size-overlap)) + 1;


step = 2; % Step between beginning of the real part of the next number
a=1;


%for m=1:(length(fftdata)/2-1)
%    audio(m) = fftdata(2*m) + j * fftdata(2*m + 1);        
%end



for n=1:time_blocks
    for m=1:freq_bins-1
        stft_image(m,n) = fftdata(a) + j * fftdata(a+1);        
        a=a+step;
    end
end
% stft_image = abs(stft_image);
stft_image = abs(stft_image);
imagesc(20*log(stft_image));
disp('Done!');

return


stft_image = abs(real(stft_image));
image(20*log(stft_image));


return
close all
clc
a=1;
for n=1:128
    for m=1:8
        hex_num = char(stft(a));
        stft_image(m,n) = hex2num(hex_num(3:end));
        hex_num = char(stft(a+1));
        stft_image(m,n) = stft_image(m,n) + j *hex2num(hex_num(3:end)); % complex part
        
        a=a+2;
    end
end
imagesc(20*log(abs(stft_image)));
disp('Done!');