%% Short test file for wavedec4 and dwt4
%
% T H   2021
clear all
close all

load('volconfig.mat'); % For easier volshow usage

% Dimensions
nx = 128;
ny = 128;
nz = 128;
nt = 64;

% Wavelet type (can also specify different wavelet for each direction)
wtype = 'db2';

% obj = randn(nx,ny,nz,nt);
% addpath('G:\N-term approx error')
obj = single(shrinkingBall(nx,ny,nz,nt));

%% Test dwt4 and idwt4

w1l = dwt4(obj,wtype);

recn1l = idwt4(w1l);

err1l = norm(obj(:)-recn1l(:))/norm(obj(:))

%% Look at stuff
k = [1,1,1,1];
wconf = 'LLLL'; wconf(k==2) = 'H';
coeff = w1l.dec(k);
figure; volshow(coeff{1}(:,:,:,1),config); title([wconf, ' wavelet coefficients'])

%% Test wavedec4 and waverec4
tic;
w = wavedec4(obj,5,wtype);

recn = waverec4(w);

err = norm(obj(:)-recn(:))/norm(obj(:))
time = toc