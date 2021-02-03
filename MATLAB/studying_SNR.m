% Studying SNR
clear
clc
close all

load('tdgf_2017_201.mat')
tdgf_c = hilbert(tdgf);

%%
Fs = 200
Ts = 1/Fs;

t = -30+Ts:Ts:30-Ts;

figure(1)
plot(t, tdgf(4000,:))
xlim([-7 0])

figure(2)
plot(t, abs(tdgf_c(4000,:)))
xlim([-7 0])