
close all
clear all
fs=8000;
secs =15;


t = 0:1/fs:2;   
trig = 5*[zeros(30,1);ones(1000,1);zeros(30,1)].';%% Trigger for oscilloscope
trig = trig - ones(1,length(trig));
%y = [trig,0.2*chirp(t,100,1,fs/2)]; 


y = randn([secs*fs 2]);
%figure;scatter(y(:,1),y(:,2));
y = [0.2,0.28; 0.4, 0.36]*y';
y =y';
%figure;scatter(y(:,1),y(:,2));
soundsc(y,fs)

                    