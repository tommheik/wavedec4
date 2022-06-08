function X = idwt4(wt,varargin)
%IDWT4 Single-level inverse discrete 4-D wavelet transform.
%   IDWT4 performs a single-level 4-D wavelet reconstruction
%   starting from a single-level 4-D wavelet decomposition.
%
%   X = IDWT4(WT) computes the single-level reconstructed 4-D array
%   X based on 4-D wavelet decomposition contained in the structure
%   WT which contains the following fields:
%     sizeINI: contains the size of the 4-D array X.
%     mode:    contains the name of the wavelet transform extension mode.
%     filters: is a structure with 4 fields LoD, HiD, LoR, HiR which
%              contain the filters used for DWT.
%     dec:     is a 2x2x2x2 cell array containing the coefficients 
%              of the decomposition.
%              dec{i,j,k,l} , i,j,k,l = 1 or 2 contains the coefficients
%              obtained by low-pass filtering (for i or j or k or l = 1)
%              or high-pass filtering (for i or j or k or l = 2).              
%
%   C = IDWT4(WT,TYPE) allows to compute the single-level reconstructed 
%   component based on the 4-D wavelet decomposition. 
%   The valid values for TYPE are:
%       - A group of 4 chars 'yxzt', one per direction, with 'y','x','z' and 't'
%         in the set {'a','d','l','h'} or in the corresponding upper case  
%         set {'A','D','L','H'}), where 'A' (or 'L') stands for low pass 
%         filter and 'D' (or 'H') stands for high pass filter.
%       - The char 'd' (or 'h' or 'D' or 'H') gives directly the sum of 
%         all the components different from the low pass one.
%
%   Examples:
%       X  = reshape(1:256,4,4,4,4);
%       wt = dwt4(X,'db1');
%       XR = idwt4(wt);
%       A  = idwt4(wt,'aaa');
%       D  = idwt4(wt,'d');
%       ADA  = idwt4(wt,'ada');
%
%   See also dwt4, wavedec4, waverec4.
%
%   Original 3D implementation by
%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Dec-2008.
%   Copyright 1995-2018 The MathWorks, Inc.
%
%   Extended for 4D by
%   T H   2021
%   University of Helsinki, Dept. of Mathematics and Statistics
%   Last edited: 12.5.2021

% Check arguments.
if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

narginchk(1,2)
s   = wt.sizeINI;
dec = wt.dec;
Lo  = wt.filters.LoR;
Hi  = wt.filters.HiR;
dwtEXTM = wt.mode;
perFLAG = isequal(dwtEXTM,'per');

lf = zeros(1,4);
for k = 1:4 , lf(k) = length(Lo{k}); end

kk = 1;
while kk<=length(varargin)
    if ischar(varargin{kk})
        word = varargin{kk};
        switch word
            case {'a','l','A','L','1'}
                                     dec{1,1,1,2}(:) = 0; 
                dec{1,1,2,1}(:) = 0; dec{1,1,2,2}(:) = 0;
                dec{1,2,1,1}(:) = 0; dec{1,2,1,2}(:) = 0; 
                dec{1,2,2,1}(:) = 0; dec{1,2,2,2}(:) = 0;
                dec{2,1,1,1}(:) = 0; dec{2,1,1,2}(:) = 0;
                dec{2,1,2,1}(:) = 0; dec{2,1,2,2}(:) = 0;
                dec{2,2,1,1}(:) = 0; dec{2,2,1,2}(:) = 0; 
                dec{2,2,2,1}(:) = 0; dec{2,2,2,2}(:) = 0;
                
            case {'d','h','D','H','0'}
                dec{1,1,1,1}(:) = 0; 
                
            otherwise
                if length(word)==4
                  num = ones(1,4);
                  for k = 1:4
                      switch word(k)
                          case {'a','l','A','L','1'}
                              num(k) = 1;
                          case {'d','h','D','H','0'} 
                              num(k) = 2;
                          otherwise
                              num(k) = -1;
                      end
                  end % for k
                  
                  for n = 1:2, for j = 1:2, for k = 1:2, for l = 1:2
                    if ~isequal([n,j,k,l],num)
                        dec{n,j,k,l}(:) = 0;
                    end
                  %  l    k    j    n
                  end, end, end, end
                else % word is not of length 4
                    error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
                end
        end
        kk = kk+1; 
    else
        s = varargin{kk};
        kk = kk+1;
    end
end

% Reconstruction. Convolve in the opposite direction:
% time steps -> slices -> rows -> columns

% Time steps
dim = 4;
U = cell(2,2,2);
for i = 1:2    
    for j = 1:2
        for k = 1:2
            U{i,j,k} = wrec1D(dec{i,j,k,1},Lo{4},dim,perFLAG,s) + ...
                 wrec1D(dec{i,j,k,2},Hi{4},dim,perFLAG,s);
        end % k
    end % j
end % i

% Slices (along z-axis)
dim = 3;
dec = cell(2,2,2); % Re-use existing variables for memory reasons
for i = 1:2    
    for j = 1:2
        dec{i,j} = wrec1D(U{i,j,1},Lo{3},dim,perFLAG,s) + ...
                 wrec1D(U{i,j,2},Hi{3},dim,perFLAG,s);
    end % j
end % i

% Rows (horizontally or along x-axis)
dim = 2;
U = cell(1,2); % Re-use existing variables for memory reasons
for i = 1:2
    U{i} = wrec1D(dec{i,1},Lo{2},dim,perFLAG,s) + ...
        wrec1D(dec{i,2},Hi{2},dim,perFLAG,s);
end

% Last reconstruction. Convolve columns (vertically or along y-axis)
dim = 1;
X = wrec1D(U{1},Lo{1},dim,perFLAG,s) + wrec1D(U{2},Hi{1},dim,perFLAG,s);
end
%-----------------------------------------------------------------------%
function Z = wrec1D(X,F,dim,perFLAG,s)

if isa(X,'single'), F = single(F); end
F = F(:); % Case to column vector

lx = size(X,dim);
lf = length(F);
nb = fix(lf/2-1);

% Permute F to "dim-dimensional" vector
switch dim
    case 1
        % Do nothing because filter is a column vector
    case 2
        F = F';
    case 3
        F = reshape(F,1,1,[]);
    case 4
        F = reshape(F,1,1,1,[]);
end

if perFLAG % Additional extension if 'per'iodic extension was used
    idxAdd = 1:nb;
    if nb>lx
        idxAdd = rem(idxAdd,lx);
        idxAdd(idxAdd==0) = lx;
    end
    idxAdd = [1:lx,idxAdd];
    lx = length(idxAdd);
    % Extend using cell array
    I = cell(1,4); I(:) = {':'};
    I{dim} = idxAdd;
end
% Odd indicies of the upsampled array
J = cell(1,4); J(:) = {':'};
J{dim} = 1:2:2*lx-1;

sZ = size(X);
if length(sZ)<4 , sZ(length(sZ)+1:4) = 1; end
sZ(dim) = 2*lx-1;

Z = zeros(sZ,class(X));
if perFLAG
    Z(J{:}) = X(I{:}); % Place elements of extended X into odd indices of Z
else
    Z(J{:}) = X; % Place elements of X into odd indices of Z
end
X = convn(Z,F);

lx = size(X,dim);
first  = floor((lx-s(dim))/2);
last  = ceil((lx-s(dim))/2);

I = cell(1,4); I(:) = {':'};

I{dim} = first + 1: lx - last;
Z  = X(I{:});
end
%-----------------------------------------------------------------------%
