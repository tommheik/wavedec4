function w = Wvec(ws)
%WVEC changes 4D wavelet structure into one long vector
%   INPUT
%   ws is the structure returned by Wavedec4.
%   OUTPUT
%   w is a long column vector.
%
%   T H   2021
A = ws.dec(1);
A = A{1}(:); % Approximation coefficients
cl = class(A);

wsz = prod(ws.sizes,2);
w = [A; zeros(15*sum(wsz),1,cl)]; % Initialize the vector
clear A
ws.dec{1} = {};
ind = wsz(1);

for l = 1:ws.level
    for k = 1:15
        di = 1+15*(l-1)+k; % dec. ind.
        coeff = ws.dec(di);
        ws.dec{di} = {};
        w(ind+1:ind+wsz(l)) = coeff{1}(:);
        ind = ind+wsz(l);
    end % k
end % l
end

