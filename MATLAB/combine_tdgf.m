function [three_year] = combine_tdgf()
%COMBINE_TDGF - takes 3 mat files for 3 seperate years and combines them
%   into a single 2D array [hours x delay]

load('tdgf_2016_201.mat');
tdgf16 = tdgf;
load('tdgf_2017_201.mat');
tdgf17 = tdgf;
load('tdgf_2018_201.mat');
tdgf18 = tdgf;
clear tdgf


three_year = vertcat(tdgf16,tdgf17,tdgf18);


end

