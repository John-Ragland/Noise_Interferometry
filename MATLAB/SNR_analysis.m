% SNR Analysis
clear
close all
clc

tdgf = combine_tdgf();

%% Plot Full Year Average

Fs = 200;
Ts = 1/Fs;

t = -30+Ts:Ts:30-Ts;

year_17 = mean(tdgf(8784:17544,:),1);

year_c = hilbert(year_17);

plot(t, abs(year_c)/max(abs(year_c)), 'linewidth',1.2)
xlim([-10,10])
ax = gca;
ax.FontSize = 18; 
title('2017 Yealong NCCF Average')
xlabel('Delay Time \tau (s)')
ylabel('Normalized Amplitude')
grid()

%% Calculate SNR for yearlong

% Look within -3.5 to -2.5 seconds
[m,bound1] = min(abs(t+2.5))
[m,bound2] = min(abs(t+3.5))

% signal_amp = max(year_17(bound2:bound1))

signal_amp = [0.8873, 1, 0.3786, 0.08623, 0.3303, 0.7609, 0.3326, 0.08309]
%Harded to look for noise between +-1 seconds
[m,idx1] = min(abs(t+1))
[m,idx2] = min(abs(t-1))

year_mag = abs(year_c);
noise_std = std(year_mag(idx1:idx2)/max(abs(year_c)))

SNR = 20*log10(signal_amp/noise_std)
%% Plot surface defined by specified windows
clear
clc

tdgf = combine_tdgf();
close all

Fs = 200;
Ts = 1/Fs;

t = -30+Ts:Ts:30-Ts;
peak_win = [5535, 5635;
    5303, 5503;
    4969, 5169;
    4622, 4782;
    6367, 6467;
    6498, 6698;
    6833, 7033;
    7242, 7402];

peak = 1;

date = 1:length(tdgf(:,1));
date = datetime(2016,1,1,0,0,0) + hours(50) + hours(date);

% Set all NAN to zero
tdgf(isnan(tdgf)) = 0;

% Create Complex Valued tdgf
% tdgf_c = hilbert(tdgf);
figure(1)
peak_name = ["A1","B1","C1","D1","A2","B2","C2","D2"];

for k = 1:8
    subplot(2,4,k)
    % surf(t(peak_win(k,1):peak_win(k,2)),date',abs(tdgf_c(:,peak_win(k,1):peak_win(k,2))));
    % shading interp
    date_st = datestr(date);
    

    
    imagesc(tdgf(:,peak_win(k,1):peak_win(k,2)));
   
    xt = get(gca, 'XTick');                         % Original 'XTick' Values
    xtlbl = linspace(t(peak_win(k,1)), t(peak_win(k,2)), numel(xt));      % New 'XTickLabel' Vector
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl);       % Label Ticks
    set(gca,'YDir','normal')
    title(peak_name(k))
    ylabel('Hours since 1/1/16 0:0:0')
    colorbar
    xlabel('Delay Time (s)')
end

%% Create SNR Plot for Each Peak
clc
clear

fprintf('\nloading tdgf...\n');
tdgf = combine_tdgf();
tdgf(isnan(tdgf))=0;
fprintf('creating complex signal...\n');
tdgf_c = hilbert(tdgf);
tdgf_mag = abs(tdgf_c);
%%
close all
Fs = 200;
Ts = 1/Fs;
t = -30+Ts:Ts:30-Ts;

[m,idx1] = min(abs(t+1.5));
[m,idx2] = min(abs(t-1.5));


peak_idx = [5585, 5403, 5069, 4702, 6417, 6598, 6933, 7322];

for k = 1:length(tdgf)
    noise = std(tdgf_mag(k,idx1:idx2));
    for n = 1:8
        amplitude(n) = tdgf(k,peak_idx(n));
        SNR(n,k) = 20*log10(amplitude(n)/noise);
    end
end

date = 1:length(tdgf(:,1));
date = datetime(2016,1,1,0,0,0) + hours(50) + hours(date);
peak_name = ["A1","B1","C1","D1","A2","B2","C2","D2"];
figure(1)
for n = 1:8
    subplot(2,4,n)
    plot(date,real(SNR(n,:)))
    ylabel('SNR (dB)')
    title(peak_name(n));

end