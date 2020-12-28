% Phase Measurement of Peaks
clear
clc
close all

load('tdgf_2017.mat');

data = tdgf_2017(1,:);
t = -30 + 1/200:1/200:30- 1/200;

idx = [5560, 5600; 5300,5500; 4900, 5200; 4680, 4720; ...
        6400, 6440; 6500, 6700; 6850, 7050; 7280, 7320];

time_instance = -1.5
[m,i] = min(abs(t-time_instance))

plot(t,data)
%% Calculate Phase for Each Peak

data_c = hilbert(data);


Data_c  =fft(data_c);
data_mag = abs(Data_c);
data_ang = angle(Data_c);

peak = 3;


% Method 1 (brute force)
data_shift = zeros(360,length(data));
for k = 1:360
    ang_new = data_ang + deg2rad(k);
    data_shift(k,:) = real(ifft(data_mag.*exp(j*ang_new)));
end
% integrate for main lag peak
integral = sum(data_shift(:,idx(peak,1):idx(peak,2)),2);

[m,phase1] = min(integral);
wrapTo360(phase1)

% Method 2 (from paper)
data_90 = real(ifft(data_mag.*exp(j*(data_ang + pi/2))));
a1 = sum(data(idx(peak,1):idx(peak,2)));
a2 = sum(data_90(idx(peak,1):idx(peak,2)));

phase2 = rad2deg(atan2(a2,a1)+pi);
wrapTo360(phase2)
figure(1)
plot(integral)

figure(2)
plot(data_shift(90,idx(peak,1):idx(peak,2)))


%% Calculate SNR for Each Peak

for peak = 1:8
    peak
    for k = 1:8658
        data_c = hilbert(tdgf_2017(k,:));

        noise_bounds = [5700, 6300]; % -1.5 and 1.5 seconds
        data_mag = abs(data_c);

        noise = std(data_mag(noise_bounds(1):noise_bounds(2)));

        SNR(peak,k) = 20*log10(max(data_mag(idx(peak,1):idx(peak,2))/noise));
    end
end

%%

plot(SNR','linewidth',1.5)
legend('Peak 1','Peak 2','Peak 3','Peak 4','Peak 5','Peak 6','Peak 7','Peak 8')