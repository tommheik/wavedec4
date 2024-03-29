function X = waverec4(wdec,varargin)
%WAVEREC4 Multilevel 4-D wavelet reconstruction.
%   WAVEREC4 performs a multilevel 4-D wavelet reconstruction
%   starting from a multilevel 4-D wavelet decomposition.
%
%   X = WAVEREC4(WDEC) reconstructs the 4-D array X based
%   on the multilevel wavelet decomposition structure WDEC.
%
%   In addition, you can use WAVEREC4 to simply extract coefficients 
%   from a 4-D wavelet decomposition (see below).
%
%   WDEC is a structure with the following fields:
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
%              'HLLL','LHLL','HHLL','LLHL','HLHL','LHHL','HHHL','LLLH',
%              'HLLH','LHLH','HHLH','LLHH','HLHH','LHHH','HHHH'. This is
%              the order in which the separable filtering operations are
%              applied in the order "columns - rows - slices - time steps"
%              NOTE: This order is different from waverec2 and waverec3
%              which use the rows -> columns convention by convolving first
%              the second dimension and then the first.
%     sizes:   contains the successive sizes of the decomposition
%              components.
%
%   C = WAVEREC4(WDEC,TYPE,N) reconstructs the multilevel
%   components at level N of a 4-D wavelet decomposition. N must be
%   a positive integer less or equal to the level of the decomposition.
%   The valid values for TYPE are:
%       - A group of 4 chars 'xyzt', one per direction, with 'x','y','z'
%         and 't' in the set {'a','d','l','h'} or in the corresponding
%         upper case  set {'A','D','L','H'}), where 'A' (or 'L') stands for
%         low pass filter and 'D' (or 'H') stands for high pass filter.
%       - The char 'd' (or 'h' or 'D' or 'H') gives directly the sum of 
%         all the components different from the low pass one.
%       - The char 'a' (or 'l' or 'A' or 'L') gives the low pass 
%         component (the approximation at level N).
%
%   For extraction purpose, the valid values for TYPE are the same as
%   above prefixed by 'c' or 'C'.
%
%   X = WAVEREC4(WDEC,'a',0) or X = WAVEREC4(WDEC,'ca',0) is equivalent
%   to X = WAVEREC4(WDEC).
%
%	C = WAVEREC4(WDEC,TYPE) is equivalent to C = WAVEREC3(WDEC,TYPE,N)
%   with N equal to the level of the decomposition.
%
%   % Example:
%       M = magic(8);
%       X = repmat(M,[1 1 8 8]);
%       wd = wavedec4(X,2,'db2','mode','per');
%       XR = waverec4(wd);
%       err1 = max(abs(X(:)-XR(:)))
%       A = waverec4(wd,'aaaa');
%       CAA = waverec4(wd,'caa');
%       D = waverec4(wd,'d');
%       err2 = max(abs(X(:)-A(:)-D(:)))
%       A1 = waverec4(wd,'aaaa',1);
%       D1 = waverec4(wd,'d',1);
%       err3 = max(abs(X(:)-A1(:)-D1(:)))
%       DDDD = waverec4(wd,'dddd',2);
%       disp(DDDD(:,:,1))
%
%   See also idwt4, wavedec4, waveinfo.
%
%   Original 3D implementation by
%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Dec-2008.
%   Copyright 1995-2018 The MathWorks, Inc.
%
%   Extended for 4D by
%   T H   2021
%   University of Helsinki, Dept. of Mathematics and Statistics
%   Last edited: 11.5.2021

% Check arguments.
if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

nbIn = length(varargin);

% Initialization.
cfs   = wdec.dec;
sizes = wdec.sizes;
level = wdec.level;
cfsFLAG = false;
nbREC_UP = level;

if nbIn>0
    errorFLAG = false;
    type = varargin{1};
    nbChar = length(type);
    % Check if first letter is 'C'
    if nbChar==2 || nbChar==3 || nbChar==5
        if isequal(upper(type(1)),'C')
            type(1) = []; nbChar = nbChar-1; cfsFLAG = true;
        else
            errorFLAG = true;
        end
    end
    % 'type' should contain one letter or four letters
    num = ones(1,4); % Decode information into num
    switch nbChar
        case 1 % 'type' is just one letter
            switch type % Check if lowpass or highpass
                case {'a','l','A','L','1'}
                    num = Inf;
                case {'d','h','D','H','0'}
                    num = -num;
                otherwise
                    errorFLAG = true;
            end
        case 4 % type is specified for each dimension
            for k = 1:4 % Go through each character
                switch type(k)
                    case {'a','l','A','L','1'}
                        % do nothing (lowpass ~ 1)
                    case {'d','h','D','H','0'}
                        num(k) = 0; % (highpass ~ 0)
                    otherwise , errorFLAG = true;
                end
            end                       
        otherwise , errorFLAG = true;
    end
    if isequal(num,[1 1 1 1])
        % if type = 'LLLL' (or similar)
        num = Inf;
    end
    
    if errorFLAG
        error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
    end
    
    if nbIn>1
        levREC = varargin{2};
        OKval = isnumeric(levREC) && isequal(levREC,fix(levREC)) && ...
            levREC<=level && (levREC>0 || ...
            (levREC==0 && (isequal(num,[1 1 1 1]) || isinf(num))) );
        if ~OKval
                error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
        end
    else
        levREC = level;
    end
    if isinf(num) % type = 'L' or 'LLLL'
        First_toDEL = 2+15*(level-levREC);
        for k = First_toDEL:(15*level+1) , cfs{k}(:) = 0; end
        if cfsFLAG
            if isequal(level,levREC) , X = cfs{1}; return; end
            nbREC_UP = level-levREC;
        end
        
    elseif isequal(num,[-1,-1,-1,-1]) % type = 'H'
        cfs{1}(:) = 0; % Remove approximation coefficients
        for k = 1:(level-levREC)
            for j = 1:15 , cfs{15*(k-1)+j+1}(:) = 0; end
        end
        if cfsFLAG , X = cfs; end
        
    else % Specific wavelet is given (for example type = 'LHHL')
        % Choose the right coefficients from right level and by decoding
        % the information in num using subbandIdx-function
        Idx_toKeep = 15*(level-levREC) + subbandIdx(num);
        if cfsFLAG
            if (~isequal(num,[1,1,1,1]) || isequal(level,levREC))
                X = cfs{Idx_toKeep};
                return
            elseif isequal(num,[1,1,1,1])
                % This shouldn't be possible because num should be set to
                % inf in this case.
                nbREC_UP = level-levREC;
            end
        end
        % Set unnecessary coefficients to 0
        for k = 1:(15*level+1)
            if k~=Idx_toKeep , cfs{k}(:) = 0; end
        end
        
    end
end

idxBeg = 1;
% Reconstruct up to level nbREC_UP
for k=1:nbREC_UP
    idxEnd = idxBeg+15;
    wdec.dec = reshape(cfs(idxBeg:idxEnd),2,2,2,2);
    X = idwt4(wdec,sizes(k+1,:)); % One-level reconstruction
    cfs{idxEnd} = X; % Save for next level
    idxBeg = idxEnd;
end
end

%-------------------------------------------------------------------------
function idxtokeep = subbandIdx(num)
% This function determines the index 'idx' of the wavelet subband given in
% 'pseudo' binary system. For example input [1 1 0 1] ~ LLHL subband should
% be 5th subband.
% 1. For the input 1=L and 0=H, but the order goes 1->0 so we revert the
%    values to match binary [1 1 0 1] -> [0 0 1 0].
% 2. While the subbands are generated by permuting the filters from right
%    to left, once the 2x2x2x2 cell (.dec) is dropped into a column the
%    order is reversed. So we by flip [0 0 1 0] -> [0 1 0 0].
% 3. This value is then converted to actual binary [0 1 0 0] -> '0100'.
% 4. Then to integer '0100' -> 4.
% 5. Finally shifted by one to match Matlab indexing 4 -> 5.
idxtokeep = 1 + bin2dec(dec2bin(fliplr(~num))');
end


