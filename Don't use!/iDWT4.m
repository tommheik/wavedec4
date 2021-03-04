function rec = iDWT4(A, D, wname, varargin)
% 4-D Discrete Wavelet Inverse Transform
%   REC = IDWT4(A,D,WNAME) returns the inverse 4-D wavelet transform of the 
%   final-level approximation coefficients A and cell array of wavelet
%   coefficients D using the wavelet WNAME.
%   A and D are outputs of DWT4. WNAME must be some real-valued wavelet
%   compatible with the WFILTERS function.
%
%   REC = IDWT4(A,D,WNAME,orgSize) computes the final synthesis so that the
%   size of REC matches the size of the original object (orgSize). This is
%   sometimes needed if the 1st level detail coefficients were excluded
%   during analysis ([A,D] = DWT4(...,"ExcludeL1")).
%
%   This code is based on the IDUALTREE4-function, just simplified.
%
%   Tommi Heikkilä
%   Created 15.5.2020
%   Last edited 18.2.2021

% A should be a 4-D matrix output from DWT4
validateattributes(A,{'numeric'},{'real','nonempty','finite','ndims',4},...
    'IDWT4','A');
% D should be a cell array of length at least one (two)
validateattributes(D,{'cell'},{'nonempty'},'IDWT4','D');

% Obtain the level of the transform
level = length(D);

% Use single precision if needed
useDouble = 1;
if isa(A, 'single')
    useDouble = 0;
end

% If original object size is not given, it is set to empty array.
if nargin < 4
    orgSize = [];
else
    orgSize = varargin{1};
end

% Check for 'adjoint' or "adjoint"
validopts = ["adjoint","inverse"];
defaultopt = "inverse";
[useAdj, ~] = ...
    wavelet.internal.getmutexclopt(validopts,defaultopt,varargin);
useAdj = strcmpi(useAdj,"adjoint"); % Use a simple boolean value

% Obtain filters
% [LoD,HiD,Lo,Hi] = wfilters(wname);
load(char("nearsym5_7"),'LoR','HiR','LoD','HiD');
Lo = LoR; Hi = HiR;

if useAdj % Use adjoint in place of inverse
    % Reconstruction filters need to be replaced with time-reversed
    % decomposition filters
    Lo = fliplr(LoD); 
    Hi = fliplr(HiD);
end

% Switch to single if needed
if ~useDouble
    Lo = single(Lo);
    Hi = single(Hi);
end

% Debug cleanup
clear HiD HiR LoD LoR

while level > 1             % Go through levels in decreasing order
    if ~isempty(D{level-1}) % Get the desired size from lower level
        syh = size(D{level-1});
        prev_level_size = syh(1:4);
    else
        syh = size(D{level});
        prev_level_size = syh(1:4).*2;
    end
    % Obtain scaling coefficients through synthesis
    A = synthesis(A,D{level},Lo,Hi,prev_level_size);
    D{level} = [];   % Clear used detail coefficients to save memory.
    level = level-1;
end

% Special case where 1st level details were excluded
if level == 1 && isempty(D{1})
    A = synthNoHighpass(A,Lo,orgSize);
end
% Synthesis from 1st level detail coefficients
if level == 1 && ~isempty(D{1})
    A = synthesis(A,D{1},Lo,Hi,[]);
end

% Return the recomposed 4-D object
rec = A;
end

%------------------------------------------------------------------------
function y = columnFilter(x,h)
% Filter the columns of x with h.
x = squeeze(x); % Get rid of nonessential dimensions
[r,c] = size(x);
y = zeros([2*r,c],class(x)); % Upsample rows by 2
y(2:2:end,:) = x; 
clear x; % Save memory

% This determines the symmetric extension of the matrix
L = length(h);
M = fix(L/2);

y = wextend('ar','sym',y,M);
if rem(L,2) == 0 % Even length filters cause trouble
    y(1,:) = []; % Remove first row
end
y = conv2(y,h(:),'valid');
end

%-------------------------------------------------------------------------
function x = synthesis(A,D,Lo,Hi,prev_level_size)
% From this we will only return the scaling coefficients at the next finer
% resolution level
Asize = size(A);
if isa(A,'single') || isa(D,'single')
    A = single(A);
    D = single(D);
    x = zeros(2*Asize,'single');
else
    x = zeros(2*Asize);
end
sy = 2*Asize;

s3a = uint8(1:Asize(3));
s4a = uint8(1:Asize(4));
s3b = Asize(3) + s3a;
s4b = Asize(4) + s4a;

%%%%%
% NOTE: THIS IMPLEMENTATION IS PROBABLY NOT THE MOST MEMORY-EFFICIENT!!!
%%%%%

% Filter dimensions in reverse order compared to analysis: X->Y->Z->T
% Filter dimensions 1 and 2
% Note that we loop through existing number of slices and time steps!
for timeidx = 1:Asize(4)
    for sliceidx = 1:Asize(3)
        % Synthesize 1st dimension
        %                   LLLL                                       HLLL
        yLLL = columnFilter(A(:,:,sliceidx,timeidx),Lo) + columnFilter(D(:,:,sliceidx,timeidx,1),Hi);
        %                   LHLL                                       HHLL             
        yHLL = columnFilter(D(:,:,sliceidx,timeidx,2),Lo) + columnFilter(D(:,:,sliceidx,timeidx,3),Hi);
        %                   LLHL                                       HLHL             
        yLHL = columnFilter(D(:,:,sliceidx,timeidx,4),Lo) + columnFilter(D(:,:,sliceidx,timeidx,5),Hi);
        %                   LHHL                                       HHHL             
        yHHL = columnFilter(D(:,:,sliceidx,timeidx,6),Lo) + columnFilter(D(:,:,sliceidx,timeidx,7),Hi);
        %                   LLLH                                       HLLH
        yLLH = columnFilter(D(:,:,sliceidx,timeidx,8),Lo) + columnFilter(D(:,:,sliceidx,timeidx,9),Hi);
        %                   LHLH                                       HHLH             
        yHLH = columnFilter(D(:,:,sliceidx,timeidx,10),Lo) + columnFilter(D(:,:,sliceidx,timeidx,11),Hi);
        %                   LLHH                                       HLHH             
        yLHH = columnFilter(D(:,:,sliceidx,timeidx,12),Lo) + columnFilter(D(:,:,sliceidx,timeidx,13),Hi);
        %                   LHHH                                       HHHH             
        yHHH = columnFilter(D(:,:,sliceidx,timeidx,14),Lo) + columnFilter(D(:,:,sliceidx,timeidx,15),Hi);
        
        % Synthesize 2nd dimension
        % Data is already upsampled in first two dimensions. For the last
        % two dimensions lowpass is stored first, then highpass.
        %
        % Each yXXX needs to be transposed before and after filtering                                        xyzt
        x(:,:,sliceidx,timeidx) = (columnFilter(yLLL.',Lo) + columnFilter(yHLL.',Hi)).'; %                   --LL
        x(:,:,Asize(3)+sliceidx,timeidx) = (columnFilter(yLHL.',Lo) + columnFilter(yHHL.',Hi)).'; %          --HL
        x(:,:,sliceidx,Asize(4)+timeidx) = (columnFilter(yLLH.',Lo) + columnFilter(yHLH.',Hi)).'; %          --LH
        x(:,:,Asize(3)+sliceidx,Asize(4)+timeidx) = (columnFilter(yLHH.',Lo) + columnFilter(yHHH.',Hi)).'; % --HH
    end
end

% Clear memory
clear A D y*

% Filter dimensions 3 and 4
% Note that we loop through doubled number of rows and columns!
for colidx = 1:sy(2)
    for rowidx = 1:sy(1)
        % Synthesize 3rd dimension
        yL = columnFilter(x(rowidx,colidx,s3a,s4a),Lo) + columnFilter(x(rowidx,colidx,s3b,s4a),Hi);
        yH = columnFilter(x(rowidx,colidx,s3b,s4a),Lo) + columnFilter(x(rowidx,colidx,s3b,s4b),Hi);
        
        % Synthesize 4th dimension
        % Each yX needs to be transposed before and after filtering 
        x(rowidx,colidx,:,:) = reshape((columnFilter(yL.',Lo) + columnFilter(yH.',Hi)).',[1,1,sy(3),sy(4)]);
    end
end
clear y*

% Now check if the size of the previous level is exactly twice the size
% of the current level. If it is exactly twice the size, the data was not
% extended at the previous level, if it is not, we have to remove the
% added row, column, and page dimensions.

if isempty(prev_level_size); return; end

% X
if prev_level_size(1) ~= sy(1)
    x = x(2:end-1,:,:,:);
end
% Y
if  prev_level_size(2) ~= sy(2)
    x = x(:,2:end-1,:,:);
end
% Z
if prev_level_size(3) ~= sy(3)
    x = x(:,:,2:end-1,:);
end
% T
if prev_level_size(4) ~= sy(4)
    x = x(:,:,:,2:end-1);
end
end

%------------------------------------------------------------------------
function y = synthNoHighpass(A,Lo,orgSize)
% We use no highpass filter here
Asize = size(A);
s1 = uint8(1:2*Asize(1));
s2 = uint8(1:2*Asize(2));
s3 = uint8(1:Asize(3));
s4 = uint8(1:Asize(4));

y = zeros(2*Asize,class(A));

% Filter dimensions in reverse order compared to analysis: X->Y->Z->T
% Filter dimensions 1 and 2
for sliceidx = s3
    for timeidx = s4
        % 1st dimension
        yl = columnFilter(A(:,:,sliceidx,timeidx),Lo);
        % 2nd dimension
        % y is tranposed before and after filtering.
        y(:,:,sliceidx,timeidx) = columnFilter(yl.',Lo).';
    end
end

% Clear memory
clear A

% Filter dimensions 3 and 4
for rowidx = s1
    for colidx = s2
        % 3rd dimension
        yl = columnFilter(squeeze(yl(rowidx,colidx,:,:)),Lo);
        % 4th dimension
        % y is transposed before and after filtering and then extended to
        % a part of 4-D array.
        y(rowidx,colidx,:,:) = reshape(columnFilter(yl.',Lo).',[1,1,2*Asize(3),2*Asize(4)]);
    end
end

% If the user specifies "excludeL1" at the input, it is possible that
% the output data size may not be correct. To correct for that, the user
% can provide the orignal data size as an input.

if ~isempty(orgSize)  % Original size is known  
    size_curr_level = size(y);
    % Dimensions are cut if needed
    if orgSize(1) ~= size_curr_level(1)
        y = y(2:end-1,:,:,:);
    end
    if  orgSize(2) ~= size_curr_level(2)
        y = y(:,2:end-1,:,:);
    end
    if orgSize(3) ~= size_curr_level(3)
        y = y(:,:,2:end-1,:);
    end
    if orgSize(4) ~= size_curr_level(4)
        y = y(:,:,:,2:end-1);
    end
end
end
