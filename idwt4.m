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
%       - A group of 4 chars 'xyzt', one per direction, with 'x','y','z' and 't'
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
%   Last edited: 3.3.2021

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
                    end
                    for n=1:2
                        for j=1:2
                            for k = 1:2
                                for l = 1:2
                                    if ~isequal([n,j,k,l],num)
                                        dec{n,j,k,l}(:) = 0;
                                    end
                                end
                            end
                        end
                    end
                else
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
% time steps -> slices -> columns -> rows

% Time steps
perm = [1,4,2,3];
U = cell(2,2,2);
for i = 1:2    
    for j = 1:2
        for k = 1:2
            U{i,j,k} = wrec1D(dec{i,j,k,1},Lo{4},perm,perFLAG,s) + ...
                 wrec1D(dec{i,j,k,2},Hi{4},perm,perFLAG,s);
        end % k
    end % j
end % i
clear dec

% Slices
perm = [1,3,2,4];
V = cell(2,2,2);
for i = 1:2    
    for j = 1:2
        V{i,j} = wrec1D(U{i,j,1},Lo{3},perm,perFLAG,s) + ...
                 wrec1D(U{i,j,2},Hi{3},perm,perFLAG,s);
    end % j
end % i
clear U

% Columns
perm = []; % Same as [1,2,3,4]
W = cell(1,2);
for i = 1:2
    W{i} = wrec1D(V{i,1},Lo{2},perm,perFLAG,s) + ...
        wrec1D(V{i,2},Hi{2},perm,perFLAG,s);
end
clear V

% Last reconstruction. Convolve rows.
perm = [2,1,3,4];
X = wrec1D(W{1},Lo{1},perm,perFLAG,s) + wrec1D(W{2},Hi{1},perm,perFLAG,s);

%-----------------------------------------------------------------------%
function X = wrec1D(X,F,perm,perFLAG,s)

if isa(X,'single'), F = single(F); end

if ~isempty(perm)
    X = permute(X,perm);
    s = s(perm);
end
if perFLAG
    lf = length(F);
    lx = size(X,2);
    nb = fix(lf/2-1);
    idxAdd = 1:nb;
    if nb>lx
        idxAdd = rem(idxAdd,lx);
        idxAdd(idxAdd==0) = lx;
    end
    X = [X X(:,idxAdd,:,:)];
end
sX = size(X);
if length(sX)<4 , sX(4) = 1; end
Z = zeros(sX(1),2*sX(2)-1,sX(3),sX(4),class(X));
Z(:,1:2:end,:,:) = X;
X = convn(Z,F);

sX = size(X,2);
F  = floor((sX-s)/2);
C  = ceil((sX-s)/2);
X  = X(:,1+F(2):end-C(2),:,:);

if ~isempty(perm) , X = ipermute(X,perm); end
%-----------------------------------------------------------------------%
