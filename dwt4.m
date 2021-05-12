function wt = dwt4(X,varargin)
%DWT4 Single-level discrete 4-D wavelet transform.
%   DWT4 performs a single-level 4-D wavelet decomposition
%   with respect to either a particular wavelet ('wname',
%   see WFILTERS for more information) or particular wavelet 
%   decomposition and reconstruction filters you specify, and 
%   using a specified DWT extension mode (see DWTMODE).
%
%   WT = DWT4(X,'wname','mode','ExtM') returns the 4-D wavelet transform
%   of the 4-D array X, 'wname' is a character vector containing the wavelet 
%   name and 'ExtM' is a character vector containing the extension mode.
%   WT = DWT4(X,'wname') uses the default extension mode: 'sym'.
%   
%   WT is a structure with the following fields:
%     sizeINI: contains the size of the 4-D array X.
%     mode:    contains the name of the wavelet transform extension mode.
%     filters: is a structure with 4 fields LoD, HiD, LoR, HiR which
%              contain the filters used for DWT.
%     dec:     is a 2x2x2x2 cell array containing the coefficients
%              of the decomposition. dec{i,j,k,l} for i,j,k,l = 1 or 2 contains
%              the coefficients obtained by lowpass filtering (for i or j
%              or k or l = 1) or highpass filtering (for i or j or k or l = 2). The
%              filtering operations are in the following order: rows,
%              columns, slices, times. For example, dec{1,2,1,1} is obtained
%              by filtering the input X along the rows with the lowpass
%              (scaling) filter, along the columns with the highpass
%              (wavelet) filter, along the slices with the
%              lowpass (scaling) filter, and along time steps with the
%              lowpass (scaling) filter.
%
%   Instead of a single wavelet, you may specify four wavelets (i.e. one
%   wavelet for each direction): WT = DWT4(X,W,...) with W =
%   {'wname1','wname2','wname3','wname4'} or W a structure with 4 fields 'w1', 'w2',
%   'w3', 'w4' containing character vectors which are the names of wavelets.
%
%   Instead of wavelets you may specify filters: 4 filters (2 for
%   decomposition and 2 for reconstruction) or 4x4 filters (one 
%   quadruplet by direction): WT = DWT4(X,WF,...)
%   Where WF must be a cell array (1x4) or (4x4) : {LoD,HiD,LoR,HiR},
%   or a structure with the four fields 'LoD', 'HiD', 'LoR', 'HiR'.
%
%   % Example:
%   %   Obtain the 4D wavelet transform of data using a wavelet name,
%   %   a cell array of filters, or a structure array.
%
%   X = reshape(1:256,4,4,4,4)
%   wt = dwt4(X,'db1')
%   [LoD,HiD,LoR,HiR] = wfilters('db2');
%   wt = dwt4(X,{LoD,HiD,LoR,HiR})
%   WS = struct('w1','db1','w2','db2','w3','db1','w4','haar');
%   wt = dwt4(X,WS,'mode','per')
%   WF = wt.filters;
%   wtBIS = dwt4(X,WF,'mode','sym')
%
%   See also dwtmode, idwt4, wavedec4, waverec4, waveinfo.

%   Original 3D implementation
%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Dec-2008.
%   Copyright 1995-2018 The MathWorks, Inc.
%   
%   Extended to 4D by 
%   T H   2021
%   University of Helsinki, Dept. of Mathematics and Statistics
%   Last edited: 11.5.2021

% Check arguments.
if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

nbIn = nargin;
narginchk(2,4);
LoD = cell(1,4); HiD = cell(1,4); LoR = cell(1,4); HiR = cell(1,4);
argStatus = true;
nextARG = 2;
if ischar(varargin{1})
    [LD,HD,LR,HR] = wfilters(varargin{1}); 
    for k = 1:4
        LoD{k} = LD; HiD{k} = HD; LoR{k} = LR; HiR{k} = HR;
    end
    
elseif isstruct(varargin{1})
    if isfield(varargin{1},'w1') && isfield(varargin{1},'w2') && ...
            isfield(varargin{1},'w3') && isfield(varargin{1},'w4')
        for k = 1:4
            [LoD{k},HiD{k},LoR{k},HiR{k}] = ...
                wfilters(varargin{1}.(['w' int2str(k)]));
        end
    elseif isfield(varargin{1},'LoD') && isfield(varargin{1},'HiD') && ...
           isfield(varargin{1},'LoR') && isfield(varargin{1},'HiR')
        for k = 1:4
            LoD{k} = varargin{1}.LoD{k}; HiD{k} = varargin{1}.HiD{k};
            LoR{k} = varargin{1}.LoR{k}; HiR{k} = varargin{1}.HiR{k};
            
        end
    else
        argStatus = false;
    end
    
elseif iscell(varargin{1})
    if ischar(varargin{1}{1})
        for k = 1:4
            [LoD{k},HiD{k},LoR{k},HiR{k}] = wfilters(varargin{1}{k});
        end
    elseif iscell(varargin{1})
        Sarg = size(varargin{1});
        if isequal(Sarg,[1 4])
            if ~iscell(varargin{1}{1})
                LoD(1:end) = varargin{1}(1); HiD(1:end) = varargin{1}(2);
                LoR(1:end) = varargin{1}(3); HiR(1:end) = varargin{1}(4);
            else
                LoD = varargin{1}{1}; HiD = varargin{1}{2};
                LoR = varargin{1}{3}; HiR = varargin{1}{4};
            end
        elseif isequal(Sarg,[3 4])
            LoD = varargin{1}(:,1)'; HiD = varargin{1}(:,2)';
            LoR = varargin{1}(:,3)'; HiR = varargin{1}(:,4)';
        else
            argStatus = false;
        end
    end
else
    argStatus = false;
end
if ~argStatus
    error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
end
sX = size(X);

% Check arguments for Extension.
dwtEXTM = 'sym';
for k = nextARG:2:nbIn-1
    switch varargin{k}
      case 'mode'  , dwtEXTM = varargin{k+1};
    end
end

dec = cell(2,2,2,2);

% Ensure that filters are column vectors so that the filtering
% operations in convn() are correct!
LoD  = cellfun(@(x)x(:),LoD,'uni',0);
LoR  = cellfun(@(x)x(:),LoR,'uni',0);
HiD  = cellfun(@(x)x(:),HiD,'uni',0);           
HiR = cellfun(@(x)x(:),HiR,'uni',0);

% Notes regarding the following: different directions of X are convolved in
% the following order: rows -> columns -> slices -> time steps.
% This is performed by permuting (reshaping) the filters because permuting
% 4D arrays is VERY inefficient.

% Convolve X across rows
dim = 1; % Dimension
% a   d  % Use bad variable names so they can be used again later
[ada,dad] = wdec1D(X,LoD{1},HiD{1},dim,dwtEXTM);

% Convolve a and d across columns
dim = 2; % Dimension
% aa  ad  % Use bad variable names so they can be used again later
[aaa,aad] = wdec1D(ada,LoD{2},HiD{2},dim,dwtEXTM);
% da  dd
[dda,ddd] = wdec1D(dad,LoD{2},HiD{2},dim,dwtEXTM);

dim = 3; % Dimension
% Convolve .. across slices (finally with correct names)
[ada,add] = wdec1D(aad,LoD{3},HiD{3},dim,dwtEXTM);
[aaa,aad] = wdec1D(aaa,LoD{3},HiD{3},dim,dwtEXTM);
[daa,dad] = wdec1D(dda,LoD{3},HiD{3},dim,dwtEXTM);
[dda,ddd] = wdec1D(ddd,LoD{3},HiD{3},dim,dwtEXTM);

dim = 4; % Dimension
% Convolve ... across time steps
[dec{1,1,1,1},dec{1,1,1,2}] = wdec1D(aaa,LoD{4},HiD{4},dim,dwtEXTM);
[dec{1,1,2,1},dec{1,1,2,2}] = wdec1D(aad,LoD{4},HiD{4},dim,dwtEXTM);
[dec{1,2,1,1},dec{1,2,1,2}] = wdec1D(ada,LoD{4},HiD{4},dim,dwtEXTM);
[dec{1,2,2,1},dec{1,2,2,2}] = wdec1D(add,LoD{4},HiD{4},dim,dwtEXTM);
[dec{2,1,1,1},dec{2,1,1,2}] = wdec1D(daa,LoD{4},HiD{4},dim,dwtEXTM);
[dec{2,1,2,1},dec{2,1,2,2}] = wdec1D(dad,LoD{4},HiD{4},dim,dwtEXTM);
[dec{2,2,1,1},dec{2,2,1,2}] = wdec1D(dda,LoD{4},HiD{4},dim,dwtEXTM);
[dec{2,2,2,1},dec{2,2,2,2}] = wdec1D(ddd,LoD{4},HiD{4},dim,dwtEXTM);

wt.sizeINI = sX;
wt.filters.LoD = LoD;
wt.filters.HiD = HiD;
wt.filters.LoR = LoR;
wt.filters.HiR = HiR;
wt.mode = dwtEXTM;
wt.dec = dec;
end

%-----------------------------------------------------------------------%
function [L,H] = wdec1D(X,Lo,Hi,dim,dwtEXTM)
% This function convolves X across dim using the vectors Lo and Hi.

if isa(X,'single'), Lo = single(Lo); Hi = single(Hi); end

lx = size(X,dim);
lf = length(Lo);

% Permute Lo and Hi to "dim-dimensional" vectors
switch dim
    case 1
        % Do nothing
    case 2
        Lo = Lo';
        Hi = Hi';
    case 3
        Lo = reshape(Lo,1,1,[]);
        Hi = reshape(Hi,1,1,[]);
    case 4
        Lo = reshape(Lo,1,1,1,[]);
        Hi = reshape(Hi,1,1,1,[]);
end

% Extend X along dim
Y = wextend4D(X,lf-1,dim,dwtEXTM);

X = convn(Y,Lo);
Z = convn(Y,Hi);

% Cutting and downsampling array (more efficient to edit index vector i
% than the whole 4D array X and Z)
lZ = size(Z,dim);
I = cell(1,4); I(:) = {':'};
% Remove extension of lf-1 from both ends
first = lf;
last = lZ-lf+1;
i = first:last;
% Remove extension caused by actual convolution
li = length(i);
first = 1+floor((li - (lx+lf-1))/2);
last = first  + (lx+lf-1) -1;
i = i(first:last);
% Downsample
first = 2;
last = length(i);
if isequal(dwtEXTM,'per')
    % Periodic extension is weird
    last = 2*ceil(lx/2);
end
i = i(first:2:last);
I{dim} = i;
% Use I to reduce the correct dimension
L = X(I{:});
H = Z(I{:});
end
%-----------------------------------------------------------------------%
function Y = wextend4D(X,a,dim,dwtEXTM)
% 4D version of wextend using correct extension for the indices
lx = size(X,dim);
% Check for allowed Extension methods
switch dwtEXTM
    case {'sym','symh','symw','sp0','per','ppd'} % Currently implemented
        % We get the indicies from the original wextend (note the limitations)
        i = wextend('ac',dwtEXTM,1:lx,a);
        % Create cell array where i is placed on dim'th place
        I = cell(1,4); I(:) = {':'};
        I{dim} = i;
        % Extend dim'th direction using i
        Y = X(I{:});
    case {'zpd','asym','asymh','asymw','sp1','spd'}
        error(['Extension method ',dwtEXTM,' not implemented!'])
        % Implementing other extension methods (zero-padding, anti-symmetric etc.)
        % would be possible by including them here (but concatenation is not
        % efficient)
    otherwise
        error('Unkown extension method.')
end
end
%-----------------------------------------------------------------------%
