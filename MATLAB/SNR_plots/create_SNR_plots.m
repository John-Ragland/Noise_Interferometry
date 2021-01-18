% Create SNR Plots for 2016, 2017, 2018

clear
clc
close all

%% 2016
load('2016_SNR_plot.mat')
% Create Datetime Variable for Y Axis
for k = 0:46
    start_hour = round(k*365/2);
    end_hour = start_hour + 365;
    mid_hour = (end_hour + start_hour)/2;
    
    year_start = datetime(2016,1,1);
    dates(k+1) = year_start + hours(mid_hour);
end

% Create avg_time variable
avg_time = 1:366;

figure1 = figure(1);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

surf(avg_time, dates, SNR)
shading interp
title('2016 SNR Plot')
ylabel('Date')
xlabel('Averaging Time (hours)')
xlim([0, 366])
view(axes1,[0 90]);

set(findall(gcf,'-property','FontSize'),'FontSize',18)
hold(axes1,'off')

%% 2017
load('2017_SNR_plot.mat')
% Create Datetime Variable for Y Axis
for k = 0:46
    start_hour = round(k*365/2);
    end_hour = start_hour + 365;
    mid_hour = (end_hour + start_hour)/2;
    
    year_start = datetime(2017,1,1);
    dates(k+1) = year_start + hours(mid_hour);
end

% Create avg_time variable
avg_time = 1:366;

figure1 = figure(1);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

surf(avg_time, dates, SNR)
shading interp
title('2017 SNR Plot')
ylabel('Date')
xlabel('Averaging Time (hours)')
xlim([0, 366])
view(axes1,[0 90]);

set(findall(gcf,'-property','FontSize'),'FontSize',18)
hold(axes1,'off')

%% 2018
load('2018_SNR_plot.mat')
% Create Datetime Variable for Y Axis
for k = 0:46
    start_hour = round(k*365/2);
    end_hour = start_hour + 365;
    mid_hour = (end_hour + start_hour)/2;
    
    year_start = datetime(2018,1,1);
    dates(k+1) = year_start + hours(mid_hour);
end

% Create avg_time variable
avg_time = 1:366;

figure1 = figure(1);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

surf(avg_time, dates, SNR)
shading interp
title('2018 SNR Plot')
ylabel('Date')
xlabel('Averaging Time (hours)')
xlim([0, 366])
view(axes1,[0 90]);

set(findall(gcf,'-property','FontSize'),'FontSize',18)
hold(axes1,'off')