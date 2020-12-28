% Analysis of Year Long TDGF Estimate
clc
clear
close all

% load tdgf.mat
load('tdgf_2016.mat');

idx = [5560, 5600; 5390,5410; 5040, 5090; 4680, 4720; ...
    6400, 6440; 6593, 6623; 6850, 7050; 7280, 7320];


% remove NaN and replace with zeros
tdgf(isnan(tdgf)) = 0;


Fs = 200; %Hz
Ts = 1/Fs;

t = -30+Ts:Ts:30-Ts;

[m, idx1] = min(abs(t+10));
[m, idx2] = min(abs(t-10));

t_crop = t(idx1:idx2);
tdgf_crop = tdgf(:,idx1:idx2);

tdgf_ds = downsample(tdgf_crop, 5);

date = 1:length(tdgf(:,1));
date = datetime(2017,1,1,0,0,0) + hours(50) + hours(date);
date = downsample(date,5);

t = -10:Ts:10;

tdgf_sym = tdgf_ds(:,2002:end) + flip(tdgf_ds(:,1:2000),2);

tdgf_complex = hilbert(tdgf_ds);

surf(abs(hilbert(tdgf(:,idx(2,1):idx(2,2)))))
shading interp
%%
surf(t,date,tdgf_ds)
datetick('y','mmm', 'keepticks')

ax = gca;
ax.FontSize = 16; 

xlabel('Propogation Time (s)')
ylabel('Month of 2017')
zlabel('Normalized Amplitude')
shading interp

%%
close all

fig1 = createfigure(t, date, abs(tdgf_complex));

exportgraphics(fig1,'yearlong_tdgf_2018.png','Resolution',500)

%%

%{

date = 1:length(tdgf(:,1));
date = datetime(2017,1,1,0,0,0) + hours(50) + hours(date);
date = downsample(date,5);
%%
close all
date = 1:length(tdgf(:,1));
date = datetime(2017,1,1,0,0,0) + hours(50) + hours(date);

[m, range_low] = min(abs(t+3.3));
[m, range_high] = min(abs(t+2.7));

[m, travel_time_idx] = max(tdgf(:,range_low:range_high), [], 2);
travel_time = t(travel_time_idx+range_low);

fig1 = figure(1)

plot(date, travel_time)

ax = gca;
ax.FontSize = 12; 

xlabel('Month')
ylabel('Propogation Time of Main Peak (s)')

% Create Moving Average

N = 721 ;
kernel = 1/N*ones(N,1);

travel_time_ma = conv(travel_time, kernel,'valid');
date_ma = date((N-1)/2:(end - (N-1)/2)-1);

figure(2)
plot(date_ma, travel_time_ma)



title('Propogation Time with 121 Point Moving Average')

ax = gca;
ax.FontSize = 12; 

xlabel('Month')
ylabel('Propogation Time of Main Peak (s)')\

%}
