% Analysis of Year Long TDGF Estimate
clc
clear
close all

% load tdgf.mat
load('tdgf_2018_201.mat');

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
peak = 2

surf(abs(hilbert(tdgf(:,idx(peak,1):idx(peak,2)))))
shading interp
%%
% Manually set invalid points
% tdgf_complex(830:1000,:) = NaN;
% Create figure
figure1 = figure;
axes1 = axes('Parent',figure1);

surf(t,date,abs(tdgf_complex)/max(abs(tdgf_complex),[],'all'))
datetick('y','mmm', 'keepticks')

axes1 = gca;
axes1.FontSize = 16; 

view(axes1,[-30.2557894736842 60]);

xlabel('Propogation Time \tau (s)')
ylabel('Month of 2018')
zlabel('Normalized Amplitude')
shading interp



%%
close all

fig1 = createfigure(t, date, abs(tdgf_complex));

exportgraphics(fig1,'yearlong_tdgf_2018.png','Resolution',500)

%% Look at Spectrogram for 2016 Spike

[s, f, t] = spectrogram(tdgf(1,:),200);

