function wdec = wavedec4(X,level,varargin)
%WAVEDEC3 Multilevel 4-D wavelet decomposition.
%   WDEC = WAVEDEC4(X,N,'wname','mode','ExtM') returns the wavelet
%   decomposition of the 4-D array X at level N, using the wavelet
%   named in string 'wname' (see WFILTERS) or particular wavelet filters
%   you specify, and using a specified DWT extension mode (see DWTMODE).
%   WDEC = WAVEDEC4(X,N,'wname') uses the default extension mode: 'sym'.
%
%   N must be a strictly positive integer (see WMAXLEV).
%
%   WDEC is the output decomposition structure, with the following fields:
%     sizeINI: contains the size of the 4-D array X.
%     level:   contains the level of the decomposition.
%     mode:    contains the name of the wavelet transform extension mode.
%     filters: is a structure with 4 fields LoD, HiD, LoR, HiR which
%              contain the filters used for DWT.
%     dec:     is an N x 1 cell array containing the coefficients of the
%              decomposition. N is equal to 15*WDEC.level+1. dec{1} contains
%              the lowpass component (approximation) at the level of the
%              decomposition. The approximation is equivalent to the
%              filtering operations 'LLLL'. dec{k+2},...,dec{k+16} with k =
%              0,15,30,...,15*(WDEC.level-1) contain the 4-D wavelet
%              coefficients. The coefficients start with the coarsest level
%              when k=0. For example, if WDEC.level=3, dec{2},...,dec{16}
%              contain the wavelet coefficients for level 3 (k=0),
%              dec{17},...,dec{31} contain the wavelet coefficients for
%              level 2 (k=15), and dec{32},...,dec{46} contain the wavelet
%              coefficients for level 1 (k=15*(WDEC.level-1)). At each
%              level, the wavelet coefficients in dec{k+2},...,dec{k+16} are
%              in the following order:
%              'HLLL','LHLL','HHLL','LLHL','HLHL','LHHL','HHHL','HLLH',
%              'LHLH','HHLH','LLHH','HLHH','LHHH','HHHH'. The strings give
%              the order in which the separable filtering operations are
%              applied from left to right. For example, suppose the order
%              is 'LHHL'. First, WAVEDEC4 applies the lowpass (scaling)
%              filter with downsampling to the rows of X. Next, WAVEDEC4
%              applies the highpass (wavelet) filter with downsampling to
%              the columns of X, then to the slices of X and then the
%              lowpass (scaling) filter to the time steps of X.
%     sizes:   contains the successive sizes of the decomposition
%              components.
%
%   % Example:
%   %   Obtain the 4D wavelet transform of a DD volume using a wavelet
%   %   name or filter pairs. Illustrate the use of the 'mode' name-value
%   %   pair.
%
%   M = magic(8);
%   X = repmat(M,[1 1 8 8]);
%   wd1 = wavedec4(X,1,'db1')
%   [LoD,HiD,LoR,HiR] = wfilters('db2');
%   wd2 = wavedec4(X,2,{LoD,HiD,LoR,HiR})
%   wd3 = wavedec4(X,2,{LoD,HiD,LoR,HiR},'mode','per')
%
%   See also dwtmode, dwt4d, waverec4, waveinfo.
%
%   Original WAVEDEC3 by
%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Dec-2008.
%   Copyright 1995-2018 The MathWorks, Inc.
%   
%   Extended to 4D by
%   T H   2021
%   Last edited: 2.3.2021

% Convert strings to char arrays
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

nbIn = nargin;
LoD = cell(1,4); HiD = cell(1,4); LoR = cell(1,4); HiR = cell(1,4);
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
            LoD = varargin{1}.LoD{k}; HiD = varargin{1}.HiD{k};
            LoR = varargin{1}.LoR{k}; HiR = varargin{1}.HiR{k};                      
        end
    else
        error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
    end
    
elseif iscell(varargin{1})
    if ischar(varargin{1}{1})
        for k = 1:4
            [LoD{k},HiD{k},LoR{k},HiR{k}] = wfilters(varargin{1}{k});
        end
    else
        LoD(1:end) = varargin{1}(1); HiD(1:end) = varargin{1}(2);
        LoR(1:end) = varargin{1}(3); HiR(1:end) = varargin{1}(4);
    end
else
    error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
end

% Check arguments for Extension.
dwtEXTM = 'sym';
for k = 2:2:nbIn-2
    switch varargin{k}
        case 'mode'  , dwtEXTM = varargin{k+1};
    end
end

% Initialization.
if isempty(X) , wdec = {}; return; end
sizes = zeros(level+1,4);
sizes(level+1,1:4) = size(X);
for k=1:level
    wdec = dwt4(X,{LoD,HiD,LoR,HiR},'mode',dwtEXTM);
    X = wdec.dec{1,1,1,1};
    if length(size(X))>2
        sizes(level+1-k,1:4) = size(X);
    else
        sizes(level+1-k,1:4) = ceil(sizes(level+2-k,1:4)/2);
    end
    wdec.dec = reshape(wdec.dec,[16,1,1,1]);
    if k>1
        cfs(1) = [];
        cfs = cat(1,wdec.dec,cfs);
    else
        cfs = wdec.dec;
    end
end
wdec.sizeINI = sizes(end,:);
wdec.level = level;
wdec.dec   = cfs;
wdec.sizes = sizes;
wdec = orderfields(wdec,{'sizeINI','level','filters','mode','dec','sizes'});

