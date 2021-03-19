clear
close all

%----------randomly simulated 3-D video data------------%
data = randi([0,255],128,128,30);

%----------Parameters------------%
% block size of (1x1x1), radius (R=1) and the number of neighbors (P=8)
R=1;
P=8;
patternMapping_u2 = getmapping(P,'u2');
nQr = 1;                   
nQc = 1; 
nQt = 1;         
rolr = 0;
colr = 0;
tolr = 0;   
%----------Parameters------------%

tic
[HIGO_xoy HIGO_xot HIGO_yot] = cal_cuboid_lbptop_vr(data, nQr, nQc, nQt, P, R, patternMapping_u2, rolr, colr, tolr);      
toc


