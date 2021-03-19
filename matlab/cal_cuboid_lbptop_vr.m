%   
%  [LBPH_xoy LBPH_xot LBPH_yot] = cal_cuboid_higo_top(Q, nQr, nQc, nQt) 
%  [LBPH_xoy LBPH_xot LBPH_yot] = cal_cuboid_higo_top(Q, nQr, nQc, nQt, P, R, lbp_mapping) 
%  [LBPH_xoy LBPH_xot LBPH_yot] = cal_cuboid_higo_top(Q, nQr, nQc, nQt, P, R, lbp_mapping, rolr, colr, tolr) 
%  Given a sequence of intensity images Q, the function returns the
%  [no_of_bins] dimensional hisogram of LBPs for three three orthogonal
%  planes, i.e. the XOY, XOT, and YOT planes.
%  
%  The input Q is arraned in a manner of [H W T], i.e., Q(:,:, t) is the t-th frame
%  of the image sequence.
% 
% nQr, nQc, and nQt stands for the numbers of divisions in row, column, and time direction 
%
% P, R, and lbp_mapping are the parameter of LBP in [0, 2*PI]. They are the
% same in the three plance but you can easily modifiy them.
%
% olr stands for overlap ratio with respect to the non-overlapping cuboid
% size.
%
% default values:
% P = 8;
% R = 1;
% lbp_mapping = getmapping(P,'u2');
% rolr = 0;
% colr = 0;
% tolr = 0;

function [LBPH_xoy LBPH_xot LBPH_yot] = cal_cuboid_lbptop(varargin)
% Version 0.99 belta
% Authors: Xiaopeng HONG, Yingyue XU and Guoying ZHAO
% Email: {xhong, yixu, gyzhao}@ee.oulu.fi
% CMV, Oulu University
%
% If you use the codes, please cite the following:
% @article{FLBPTOP_2016,
%     title    = {LBP-TOP: a Tensor Unfolding Revisit},
%     author   = {Xiaopeng Hong and Yingyue Xu and Guoying Zhao},
%     journal  = {ACCV Workshop on "Spontaneous Facial Behavior Analysis"},
%     year     = {2016}
% }
%
% Changelog
% The implementation of LBP-TOP is based on the following paper:
% Zhao G & Pietik?inen M (2007) Dynamic texture recognition using local binary patterns with an application to facial expressions. IEEE Transactions on Pattern Analysis and Machine Intelligence, 29(6):915-928.
% NOTE THAT: this code is for general purpose hence the division from sequence
% to cuboids is just an example. You can choose whatever division strategy
% you like. Moreover, the P and R in each of the orthogonal plances are not
% necessary to be the same as we implement here.
%
%ALSO NOTE THAT: we DID NOT cut the two frames in the beginning and the end
%of the sequence to maintain as many frames as we can get, hence we padd the sequence in the Temporal direction, you
%may choose whether to pad or to cut in the Temporal direction according to
%your applications.
% in the case without padding, use divisons_Arr = get_cuboid_divisons(hs-2,
% ws-2, fs-2, nQr, nQc, nQt, rolr, colr, tolr); and revising all statements of (fs+2)
% 

if nargin~=4 && nargin~=7 && nargin~=10
    error('Invalid number of arguments');
    LBPH_xoy = [];
    LBPH_xot = [];
    LBPH_yot = [];
    return;
end

Q = varargin{1};

nQr = varargin{2};
nQc = varargin{3};
nQt = varargin{4};
% default setting for the remaining parameters
%no_of_bins = 4;
% P = varargin{5};
% lbp_mapping = getmapping(P,'u2');
%
if nargin>=5
P = varargin{5}; 
R = varargin{6};  
lbp_mapping = varargin{7};  
end
% olr stands for overlap ratio with respect to the non-overlapping cuboid
% size.
%
if nargin == 10
rolr = varargin{8};
colr = varargin{9};
tolr = varargin{10};
end

no_of_bins = lbp_mapping.num;

%% TO Divide the sequence into cuboids

[hs ws fs] = size(Q);
divisons_Arr = get_cuboid_divisons(hs-2*R, ws-2*R, fs-2*R, nQr, nQc, nQt, rolr, colr, tolr);
% check the number of cuboid divisons
[coordinate_notused nQr2 nQc2 nQt2] = size(divisons_Arr);
assert(nQr2 == nQr)
assert(nQc2 == nQc)
assert(nQt2 == nQt)


%% calculate LBP & distributions in the XOY planes

I_HWT = zeros(hs, 2*R + ws *fs);
I_HWT(:, R+1:end-R) = reshape(Q, [hs, ws*fs]);

%call any meaningful implementation of lbp in 2D image...
QLBP_xoy =LBP(I_HWT, R, P, lbp_mapping,'I');
QLBP_xoy = reshape(QLBP_xoy, [hs-2*R, ws, fs]);
QLBP_xoy = QLBP_xoy(:,1+R:end-R,1+R:end-R);

%%  And ... in the YOT plane
Q_THW = permute(Q, [3 1 2]);    %T H W
I_THW =  zeros(fs, 2*R + ws *hs);
I_THW(:, 1+R:end-R) = reshape(Q_THW, [size(Q_THW,1), size(Q_THW,2)*size(Q_THW,3)]);
QLBP_yot =LBP(I_THW, R, P,lbp_mapping,'I');
QLBP_yot = reshape(QLBP_yot, [fs-2*R, hs, ws]);
QLBP_yot = permute(QLBP_yot, [2 3 1]);    %H T W

QLBP_yot = QLBP_yot(1+R:end-R,1+R:end-R,:);

%% ... And in the XOT plane
Q_TWH = permute(Q, [3 2 1]);    %T W H
I_TWH =  zeros(fs, 2*R + ws *hs);
I_TWH(:, 1+R:end-R) = reshape(Q_TWH, [size(Q_TWH,1), size(Q_TWH,2)*size(Q_TWH,3)]);
QLBP_xot = LBP(I_TWH, R, P,lbp_mapping,'I');
QLBP_xot = reshape(QLBP_xot, [fs-2*R, ws, hs]);
QLBP_xot = permute(QLBP_xot, [3 2 1]);    %H T W

QLBP_xot = QLBP_xot(1+R:end-R,1+R:end-R,:);


%% to get the histogram from three planes

LBPH_xoy = zeros(no_of_bins, nQr*nQc*nQt);
LBPH_yot = zeros(no_of_bins, nQr*nQc*nQt);
LBPH_xot = zeros(no_of_bins, nQr*nQc*nQt);

divisons_Arr = reshape(divisons_Arr, [6, numel(divisons_Arr)./6]);

for divNo = 1 : size(divisons_Arr,2)
    
    %get the current cuboid division here
    cur_div = divisons_Arr(:, divNo);
    r_beg = cur_div(1); r_end = cur_div(2);
    c_beg = cur_div(3); c_end = cur_div(4);
    t_beg = cur_div(5); t_end = cur_div(6);
    
    LBP_cuboid = QLBP_xoy(r_beg:r_end,c_beg:c_end, t_beg:t_end);
    h = hist(LBP_cuboid(:)+1, 1:no_of_bins); %check the index according to your getmapping function
    LBPH_xoy(:,divNo) = h./sum(h);
    
    LBP_cuboid = QLBP_xot(r_beg:r_end,c_beg:c_end, t_beg:t_end);
    h = hist(LBP_cuboid(:)+1, 1:no_of_bins); %check the index according to your getmapping function
    LBPH_xot(:,divNo) = h./sum(h);
    
    LBP_cuboid = QLBP_yot(r_beg:r_end,c_beg:c_end, t_beg:t_end);
    h = hist(LBP_cuboid(:)+1, 1:no_of_bins); %check the index according to your getmapping function
    LBPH_yot(:,divNo) = h./sum(h);
    
end

LBPH_xoy = LBPH_xoy(:)';
LBPH_xot = LBPH_xot(:)';
LBPH_yot = LBPH_yot(:)';

return;
end

function vec = arr2vec(mat)
    vec = mat(:);
end
function [divisons_Arr] = get_cuboid_divisons(varargin)
% [divisons_Arr] = get_cuboid_divisons(R, C, T, nr, nc, nf)
% [divisons_Arr] = get_cuboid_divisons(R, C, T, nr, nc, nf, rolr, colr, tolr)
%
% the function returns a matrix that contains the begining points and the ending points in three dimensions (row/y, column/x, and t) for each of the [nr by nc by nf] cuboids.
% size(ivisons_Arr) = [6, nr, nc, nf];
% divisons_Arr(:, m, n, l) = [r_beg(m),r_end(m), c_beg(n), c_end(n), t_beg(l), t_end(l)]
% R C and T stands for the height, width, and frame numbers of the image
% sequence to divide
%
% nr, nc, and nf stands for the numbers of division in row, column, and time direction 
%
% olr stands for overlap ratio with respect to the non-overlapping cuboid
% size.
%
% default values:
% rolr = 0;
% colr = 0;
% tolr = 0;
%% VERY IMPORTANT !!!
% I have been trying different combinations of the functions changing floating
% points to integers
% if anyone have ideas of how to improve the get_cuboid_divisons()
% function, please kindly let me know: xhong@ee.oulu.fi. Thank you in
% advance!

if nargin~=6&&nargin~=9
    error('Invalid number of arguments');
end

R = varargin{1};
C = varargin{2};
T = varargin{3};

nr = varargin{4};
nc = varargin{5};
nf = varargin{6};

rolr = 0;
colr = 0;
tolr = 0;

if nargin==9
    rolr = varargin{7};
    colr = varargin{8};
    tolr = varargin{9};
end


if R < nr
    error('Number of row blocks larger than the image dimension!');
end

if C < nc
    error('Number of column blocks larger than the image dimension!');
end

if T < nf
    error('Number of time blocks larger than the image dimension!');
end

NooverBlockW = ceil(C/nc);
NooverBlockH = ceil(R/nr);
NooverBlockT = ceil(T/nf);

% NooverBlockW = round(C/nc);
% NooverBlockH = round(R/nr);
% NooverBlockT = round(T/nf);

XoverSize = floor(NooverBlockW*colr);
YoverSize = floor(NooverBlockH*rolr);
ToverSize = floor(NooverBlockT*tolr);

BlockW = round((C+XoverSize*(nc-1))/nc);
BlockH = round((R+YoverSize*(nr-1))/nr);
BlockT = round((T+ToverSize*(nf-1))/nf);


stepX = max(NooverBlockW - XoverSize, 1);
stepY = max(NooverBlockH - YoverSize, 1);
stepT = max(NooverBlockT - ToverSize, 1);

BlockW = max(stepX, BlockW);
BlockH = max(stepY, BlockH);
BlockT = max(stepT, BlockT);
%stepT = min(stepT, BlockT);

c_beg = 1:stepX:C-BlockW+1;
if numel(c_beg) > nc
    assert(numel(c_beg) == nc+1, ['nc+1=' num2str(nc+1) ' real nc=' num2str(numel(c_beg))])
    c_beg = c_beg(1:nc);
    c_beg(nc) = c_beg(nc-1) + stepX;    %min(c_beg(nc-1) + stepX, C - BlockW+1);
    
elseif numel(c_beg) < nc
    assert(numel(c_beg) == nc-1)
    c_beg(nc) = c_beg(nc-1) + stepX;    %max(c_beg(nc-1) + stepX, C - BlockW+1);
end

r_beg = 1:stepY:R-BlockH+1;
if numel(r_beg) > nr
    assert(numel(r_beg) == nr+1, ['nr+1=' num2str(nr+1) ' real nr=' num2str(numel(r_beg))])
    r_beg = r_beg(1:nr);
    r_beg(nr) = r_beg(nr-1) + stepY;    %min(r_beg(nr-1) + stepY, R - BlockH+1);
elseif numel(r_beg) < nr
    assert(numel(r_beg) == nr-1)
    r_beg(nr) = r_beg(nr-1) + stepY;    %max(r_beg(nr-1) + stepY, R - BlockH+1);
end

t_beg = 1:stepT:T-BlockT+1;
if numel(t_beg) > nf
    assert(numel(t_beg) == nf+1, ['nf+1=' num2str(nf+1) ' real nf=' num2str(numel(t_beg))])
    t_beg = t_beg(1:nf);
    t_beg(nf) = t_beg(nf-1) + stepT;    %min(t_beg(nf-1) + stepT, T - BlockT+1);    
elseif numel(t_beg) < nf    
    assert(numel(t_beg) == nf-1)
    t_beg(nf) = t_beg(nf-1) + stepT;    %max(t_beg(nf-1) + stepT, T - BlockT+1); 
end

c_end = c_beg + BlockW - 1; c_end(end) = C;
r_end = r_beg + BlockH - 1; r_end(end) = R;
t_end = t_beg + BlockT - 1; t_end(end) = T;

assert(nr == numel(r_end), ['nr=' num2str(nr) ' real nr=' num2str(numel(r_end))]);
assert(nc == numel(c_end), ['nc=' num2str(nc) ' real nc=' num2str(numel(c_end))]);
assert(nf == numel(t_end), ['nf=' num2str(nf) ' real nf=' num2str(numel(t_end))]);
divisons_Arr = zeros(6, nr, nc, nf);
for m = 1 : nr
    
    for n = 1 : nc
        
        for l = 1 : nf
            divisons_Arr(:, m, n, l) = [r_beg(m),r_end(m), c_beg(n), c_end(n), t_beg(l), t_end(l)];
        end
    end
end
end


