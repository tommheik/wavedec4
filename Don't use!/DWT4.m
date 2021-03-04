function [A, D] = DWT4(f, level, wname, useDouble, varargin)
% 4-D Discrete Wavelet Transform
%   [A,D] = DWT(F) returns the 4-D discrete complex wavelet
%   transform of F at the maximum level, floor(log2(min(size(F)))). F is a
%   real-valued 4-D array (X-by-Y-by-Z-by-T) where all dimensions (X,Y,Z,T)
%   must be even and greater than or equal to 4. By default, Haar wavelets
%   are used but most usual wavelet filters are available.
%   A is the matrix of real-valued final-level
%   scaling (lowpass) coefficients. D is a 1-by-L cell array of wavelet
%   coefficients, where L is the level of the transform. There are 15
%   wavelet subbands in the 4-D DWT at each level. The wavelet coefficients 
%   are real-valued. 
%
%   [A,D] = DWT4(X,LEVEL) obtains the 4-D discrete wavelet transform down
%   to LEVEL. LEVEL is a positive integer greater than or equal to 2 and
%   less than or equal to floor(log2(min(size(F))).
%   
%   WNAME is the name of the wavelet to be used. Valid options are at least
%   'haar', 'db*' (for * = 1,2,...,45), symlets etc. Check WFILTERS for
%   more information.
%
%   USEDOUBLE toggles between double precision arrays (useDouble = 1) and
%   singles precision (useDouble = 0). Class of f is used by default.
%
%   "ExcludeL1" excludes the first level detail coefficients and only the
%   lowpass filter is used. If this option is used a perfect reconstruction
%   is no longer possible but both the forward and inverse transform become
%   computationally more efficient.
%
%   This code is based on the DUALTREE4-function, just simplified.
%
%   Tommi Heikkilä
%   Created 16.2.2021
%   Last edited 17.2.2021

% Ensure the input is numeric, real, and that it is four-dimensional
validateattributes(f,{'numeric'},{'real','nonempty','finite','ndims',4},...
    'DWT4','F');

% Check that all the dimensions of x are even and every dimension is
% greater than or equal to 4
origsizedata = size(f);
if any(rem(origsizedata,2)) || any(origsizedata < 4)
    error('Object dimensions are incompatible');
end
% Set up defaults for the lowpass and highpass filters
maxlev = floor(log2(min(size(f))));
if nargin >= 2
validateattributes(level,{'numeric'},...
        {'integer','scalar','<=',maxlev,'>=',1},'DWT4','LEVEL');
    params.level = level;
end

if nargin < 2 % Default to maximum level
    level = maxlev;
end

if nargin < 3 % Default to Haar wavelets
    wname = 'haar';
end

if nargin < 4 % Default to whatever class the input is
    useDouble = isa(f,'double');
end

% Check for 'ExcludeLeve1Details' or "ExcludeLevel1Details"
validopts = ["ExcludeL1","IncludeL1"];
defaultopt = "IncludeL1";
[opt, varargin] = ...
    wavelet.internal.getmutexclopt(validopts,defaultopt,varargin);

% Case f to double or single
if useDouble
    f = double(f);
else
    f = single(f);
end

% Obtain the analysis filters
% [Lo,Hi] = wfilters(wname);
load(char("nearsym5_7"),'LoD','HiD');
Lo = LoD; Hi = HiD;
if ~useDouble
    Lo = single(Lo);
    Hi = single(Hi);
end

% Allocate array for wavelet coefficients
D = cell(level,1);

% Level 1 filtering. We can omit the highest level
% details if needed
if strcmpi(opt,"ExcludeL1")
    A = analysisNoHighpass(f,Lo);
    D{1} = [];
else
    [A,D{1}] = analysisHighpass(f,Lo,Hi);
end

lev = 2;
% For levels two and up, we always keep detail coefficients
while lev <= level    
    [A,D{lev}] = analysisHighpass(A,Lo,Hi);
    lev = lev+1;
end
end

%------------------------------------------------------------------------
function y = columnFilter(x,h)
% Filter the columns of x with h. This function does not decimate the
% output.

% This determines the symmetric extension of the matrix
L = length(h);
M = fix(L/2);

x = wextend('ar','sym',x,M);
if rem(L,2) == 0 % Even length filters cause trouble
    x(1,:) = []; % Remove first row
end
y = conv2(x,h(:),'valid');
y = y(1:2:end,:); % Downsample by 2
end

%-------------------------------------------------------------------------
function A = analysisNoHighpass(x,Lo)
% This function is called if the user specified "excludeL1"
sx = size(x);

% Filter dimensions 4 and 3
for rowidx = 1:sx(1)
   for colidx = 1:sx(2)
       % Select one 2D array and permute it so that the 4th dimension
       % forms the columns.
       %                                              t z x y
       y = columnFilter(permute(x(rowidx,colidx,:,:),[4,3,1,2]),Lo);
       % Transpose and filter dimension 3.
       x(rowidx,colidx,:,:) = reshape(columnFilter(y.',Lo),[1,1,sx(3),sx(4)]);
   end
end

% Filter dimensions 2 and 1
for sliceidx = 1:sx(3)
   for timeidx = 1:sx(4)
       % Select one 2D array and transpose it so that the 2nd dimension
       % forms the columns.
       y = columnFilter(x(:,:,sliceidx,timeidx).',Lo);
       % Transpose and filter dimension 1.
       x(:,:,sliceidx,timeidx) = columnFilter(y.',Lo);
   end
end

A = x;
end

%------------------------------------------------------------------------
function [A,D] = analysisHighpass(x,Lo,Hi)
% This function computes one decomposition level

sx = size(x);
hpx = ceil(sx/2); % Half point of every dimension
xtmp = zeros(2*hpx,class(x)); % xtmp might be slightly larger than x if x is of odd length

% Note this is half of the original input size
% "lowpass" indicies
s1a = uint8(1:hpx(1));
s2a = uint8(1:hpx(2));
s3a = uint8(1:hpx(3));
s4a = uint8(1:hpx(4));

% Latter half 
% "highpass indicies"
s1b = hpx(1)+s1a;
s2b = hpx(2)+s2a;
s3b = hpx(3)+s3a;
s4b = hpx(4)+s4a;

% Filter dimensions 4 and 3
% Note that we only loop through original number of rows and columns!
for rowidx = 1:sx(1)
   for colidx = 1:sx(2)
       % Select one 2D array and permute it so that the 4th dimension
       % forms the columns.
       %                                 t z x y
       y = permute(x(rowidx,colidx,:,:),[4,3,1,2]); 
       %       Lowpass              Highpass
       y = [columnFilter(y,Lo); columnFilter(y,Hi)].'; % Combine and transpose
       sy = [1,1,hpx(3),2*hpx(4)]; % Consider size as a part of 4-D array
       
       % 3rd dimension
       xtmp(rowidx,colidx,s3a,:) = reshape(columnFilter(y,Lo),sy); % Lowpass
       xtmp(rowidx,colidx,s3b,:) = reshape(columnFilter(y,Hi),sy); % Highpass
   end
end

% Save some memory
clear x

% Filter dimensions 2 and 1
% Note that we loop through possibly slighly larger number of slices and time steps!
for sliceidx = 1:2*hpx(3)
    for timeidx = 1:2*hpx(4)
        % Select one 2D array and permute it so that the 2nd dimension
        % forms the columns.
        %                                           y x z t
        y = permute(xtmp(:,:,sliceidx,timeidx),[2,1,3,4]); 
        %       Lowpass              Highpass
        y = [columnFilter(y,Lo); columnFilter(y,Hi)].'; % Combine and transpose
        sy = [1,1,hpx(1),2*hpx(2)]; % Consider size as a part of 4-D array
        
        % 1st dimension
        xtmp(s1a,:,sliceidx,timeidx) = reshape(columnFilter(y,Lo),sy); % Lowpass
        xtmp(s1b,:,sliceidx,timeidx) = reshape(columnFilter(y,Hi),sy); % Highpass
    end
end

% Note in listing the subbands the order is reversed compared to what was
% done previously, i.e. 1st dimension, then 2nd, then 3rd and then 4th.
A = xtmp(s1a,s2a,s3a,s4a);                   % LLLL
% Save the different detail coefficients
D = cat(5, xtmp(s1b,s2a,s3a,s4a),...         % HLLL
           xtmp(s1a,s2b,s3a,s4a),...         % LHLL
           xtmp(s1b,s2b,s3a,s4a),...         % HHLL
           xtmp(s1a,s2a,s3b,s4a),...         % LLHL
           xtmp(s1b,s2a,s3b,s4a),...         % HLHL
           xtmp(s1a,s2b,s3b,s4a),...         % LHHL
           xtmp(s1b,s2b,s3b,s4a),...         % HHHL
           xtmp(s1a,s2a,s3a,s4b),...         % LLLH
           xtmp(s1b,s2a,s3a,s4b),...         % HLLH
           xtmp(s1a,s2b,s3a,s4b),...         % LHLH
           xtmp(s1b,s2b,s3a,s4b),...         % HHLH
           xtmp(s1a,s2a,s3b,s4b),...         % LLHH
           xtmp(s1b,s2a,s3b,s4b),...         % HLHH
           xtmp(s1a,s2b,s3b,s4b),...         % LHHH
           xtmp(s1b,s2b,s3b,s4b));           % HHHH
end
