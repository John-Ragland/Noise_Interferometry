% create SNR Plots 2

clear
clc
close all
load('2017_SNR_plot.mat')

for k = 1:8
    subplot(2,4,k)
    imagesc(squeeze(SNRs(k,:,:)))
end