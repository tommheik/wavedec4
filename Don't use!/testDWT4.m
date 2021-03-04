% Simple test file for DWT4
%
% T H    2021
clear all
close all

x = 64;
y = 64;
z = 64;
t = 64;

level = 4;
wname = 'sym7';

vecNormSqr = @(x) norm(x(:))^2;

addpath('G:\N-term approx error\')
sb = single(shrinkingBall(x,y,z,t));

sbNorm = vecNormSqr(sb);

[A,D] = DWT4(sb,level,wname);

figure; montage(A(:,:,:,1));

coeffNorm = vecNormSqr(A) + sum(cellfun(vecNormSqr, D));

sbRecn = iDWT4(A,D,wname);

figure; montage(sbRecn(:,:,:,1));

recnNorm = vecNormSqr(sbRecn);