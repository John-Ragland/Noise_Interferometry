% create SNR Plots 2

clear
clc
close all

% Create Datetime Variable for Y Axis
for k = 0:46
    start_hour = round(k*365/2);
    end_hour = start_hour + 365;
    dates(k+1) = (end_hour + start_hour)/2;
end

avg_time = 1:366;
load('2017_SNR_plot.mat')
peak_names = ["dA", "s1b0A", "s2b1A", "s3b2A", "dB", "s1b0B", "s2b1B", "s3b2B"]

bottom = min(SNRs,[],'all')
top = max(SNRs,[],'all')
figure(1)
for k = 1:8
    subplot(2,4,k)
    imagesc(avg_time,dates, squeeze(SNRs(k,:,:)))
    set(gca,'YDir','normal')
    title(peak_names(k))
    xlabel('Average Hours')
    ylabel('Hours of 2017')
    caxis manual
    caxis([bottom top]);

end
set(findall(gcf,'-property','FontSize'),'FontSize',14)
hp8 = get(subplot(2,4,8),'Position');
hp4 = get(subplot(2,4,4),'Position');
colorbar('Position', [hp8(1)+hp8(3)+0.04  hp8(2)  0.02  hp4(2)+hp4(3)*1.5])

%%
% Average SNR Plot
figure(2)

% set all NAN to 0
nan_mask = isnan(SNRs);
SNRs(nan_mask) = 0;
for k = 1:8
    subplot(2,4,k)
    
    snr_avg = mean(squeeze(SNRs(k,:,:)),1);
    plot(snr_avg, 'linewidth',1.5)
    title(peak_names(k))
    xlabel('Average Hours')
    ylabel('SNR (dB)')
    caxis manual
    caxis([bottom top]);

end

set(findall(gcf,'-property','FontSize'),'FontSize',14)
