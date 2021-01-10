% Air Gun Analysis

load('tdgf_july-aug.mat');

tdgf(isnan(tdgf)) = 0;

tdgf_c = hilbert(tdgf);
Fs = 200;
Ts = 1/Fs
t = -30+Ts:Ts:30-Ts;
surf(abs(tdgf_c))
shading interp