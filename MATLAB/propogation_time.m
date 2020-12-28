clear
clc
close all

% Propogation Time Analysis
Fs = 200;
Ts = 1/Fs;
t = -30+Ts:Ts:30-Ts;

load('tdgf_2017.mat');

date = 1:length(tdgf(:,1));
date = datetime(2017,1,1,0,0,0) + hours(50) + hours(date);

[m, range_low] = min(abs(t+3.3));
[m, range_high] = min(abs(t+2.7));

[m, travel_time_idx] = max(tdgf(:,range_low:range_high), [], 2);
travel_time_17 = t(travel_time_idx+range_low);

load('tdgf_2016.mat')
tdgf = tdgf_2016;
[m, travel_time_idx] = max(tdgf(:,range_low:range_high), [], 2);
travel_time_16 = t(travel_time_idx+range_low);

fig1 = figure(1)

plot(date, travel_time_16)
hold on
plot(date, travel_time_17)
hold off

ax = gca;
ax.FontSize = 12; 

xlabel('Month')
ylabel('Propogation Time of Main Peak (s)')

% Create Moving Average
N = 721 ;
kernel = 1/N*ones(N,1);

travel_time_ma_17 = conv(travel_time_17, kernel,'valid');
travel_time_ma_16 = conv(travel_time_16, kernel,'valid');
date_ma = date((N-1)/2:(end - (N-1)/2)-1);

figure(2)
plot(date_ma, travel_time_ma_16,'linewidth',2)
hold on
plot(date_ma, travel_time_ma_17,'linewidth',2)
legend('2016','2017')
title('Propogation Time with 121 Point Moving Average')

ax = gca;
ax.FontSize = 12; 

xlabel('Month')
ylabel('Propogation Time of Main Peak (s)')

%% peak interpolation

clear
clc
close all

peak = 2;
% Propogation Time Analysis
Fs = 200;
Ts = 1/Fs;
t = -30+Ts:Ts:30-Ts;

load('tdgf_2017.mat');
idx = [5560, 5600; 5300,5500; 4900, 5200; 4680, 4720; ...
        6400, 6440; 6500, 6700; 6850, 7050; 7280, 7320];
    
tdgf_c = hilbert(tdgf_2017);
tdgf_mag = abs(tdgf_c);

[m, travel_time_idx] = max(tdgf_mag(:,idx(peak,1):idx(peak,2)), [], 2);

%%

i = travel_time_idx + idx(2,1) - 1;

for k = 1:8658
    alpha(k) = tdgf_mag(k,i(k)-1);
    beta(k) = tdgf_mag(k,i(k));
    gamma(k) = tdgf_mag(k,i(k)+1);
end

%%
p = 1/2*(alpha - gamma)./(alpha - 2.*beta + gamma);

bin = i + p';
prop_time = Ts*bin - 30;

%% Moving Average

% Create Moving Average
date = 1:length(tdgf_mag(:,1));
date = datetime(2017,1,1,0,0,0) + hours(50) + hours(date);

N = 1;
kernel = 1/N*ones(N,1);

travel_time_ma = conv(prop_time, kernel,'valid');
if N == 1
    date_ma = date;
else 
    date_ma = date((N-1)/2:(end - (N-1)/2)-1);
end

figure(2)
plot(date_ma, travel_time_ma)



title('Propogation Time with 121 Point Moving Average')

ax = gca;
ax.FontSize = 12; 

xlabel('Month')
ylabel('Propogation Time of Main Peak (s)')


%% Use function to create figure
clear
clc

[date, prop_16, SNR16] = ma_prop_time('tdgf_2016.mat', 1, 2);
[date, prop_17, SNR17] = ma_prop_time('tdgf_2017.mat', 1, 2);
[date, prop_18, SNR18] = ma_prop_time('tdgf_2018.mat', 1, 2);

%%
close all
fig1 = figure(1)
ax = gca;
ax.FontSize = 14; 
figure(1)
hold on
plot(date, prop_16, 'linewidth',2)
plot(date, prop_17, 'linewidth',2)
plot(date, prop_18, 'linewidth',2)
hold off
hold off

% Add y tick formatting
datetick('x','mmm', 'keepticks')

legend('2016','2017','2018')
grid()
ylabel('Propagation Time (s)')
xlabel('Month')

% exportgraphics(fig1,'prop_time_3yr_1.png','Resolution',500)

fig2 = figure(2)
plot(date, SNR16, 'linewidth',1.5)
xlabel('Month')
datetick('x','mmm', 'keepticks')
ylabel('SNR (dB)')
grid()
title('SNR for 2016')

exportgraphics(fig2,'SNR16.png','Resolution',500)

fig3 = figure(3)
plot(date, SNR17, 'linewidth',1.5)
xlabel('Month')
datetick('x','mmm', 'keepticks')
ylabel('SNR (dB)')
grid()
title('SNR for 2017')
exportgraphics(fig3,'SNR17.png','Resolution',500)

fig4 = figure(4)
plot(date, SNR18, 'linewidth',1.5)
xlabel('Month')
datetick('x','mmm', 'keepticks')
ylabel('SNR (dB)')
grid()
title('SNR for 2018')
exportgraphics(fig4,'SNR18.png','Resolution',500)
%%

close all
plot(prop_16)
%%
plot(prop_16)