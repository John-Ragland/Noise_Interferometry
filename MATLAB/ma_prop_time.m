function [date_vector, prop_time_ma, SNR] = ma_prop_time(filename,ma_size, peak_number)
%MOVING_AVERAGE_PROP_TIME Finds movings average of the propogation time
%between two hydrophones using noise interferometry experiment.
%   - reads noise_interferometry experiment saved in filename
%   - uses parabolic interpolation to find time of peak for each hour in
%       year. For peak number peak_number (convention to specified later)
%   - computes moving average of size ma_size to propogation time vector
%
%Returns:
%   - date_vector: vector with time and date of each sample
%   - prop_time_ma: estimate from tdgf of propogation time
%   - SNR: Signal to noise raitio of peak for time slot

Fs = 200;
Ts = 1/Fs;

%idx = [5560, 5600; 5300,5500; 4900, 5200; 4680, 4720; ...
%    6400, 6440; 6500, 6700; 6850, 7050; 7280, 7320];

noise_range = [5700, 6300];

idx = [5560, 5600; 5400,5405; 5040, 5090; 4680, 4720; ...
    6400, 6440; 6593, 6623; 6850, 7050; 7280, 7320];


load(filename);
% remove NaN and replace with zeros
tdgf(isnan(tdgf)) = 0;
tdgf_c = hilbert(tdgf);
tdgf_mag = abs(tdgf_c);

[m, travel_time_idx] = max(tdgf_mag(:,idx(peak_number,1):idx(peak_number,2)), [], 2);

i = travel_time_idx + idx(2,1) - 1;
alpha = zeros(8658,1);
beta = zeros(8658,1);
gamma = zeros(8658,1);
for k = 1:8658
    alpha(k) = tdgf_mag(k,i(k)-1);
    beta(k) = tdgf_mag(k,i(k));
    gamma(k) = tdgf_mag(k,i(k)+1);
end

p = 1/2*(alpha - gamma)./(alpha - 2.*beta + gamma);

bin = i + p;
prop_time = Ts*bin - 30;

% Moving Average

% Create Moving Average
date = 1:length(tdgf_mag(:,1));
date = datetime(2017,1,1,0,0,0) + hours(50) + hours(date);

N = ma_size;
kernel = 1/N*ones(N,1);

prop_time_ma = conv(prop_time, kernel,'valid');
if N == 1
    date_vector = date;
else 
    date_vector = date((N-1)/2:(end - (N-1)/2)-1);
end
    
    
noise_std = std(tdgf(:,noise_range(1):noise_range(2)),0,2);

for k = 1:8658
    peak_amp(k) = tdgf_mag(k,i(k));
end

SNR = 20*log10(peak_amp./noise_std');

end

