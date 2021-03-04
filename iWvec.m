function ws = iWvec(w,ws)
%IWVEC changes long wavelet vector into one 4D wavelet structure
%   INPUT
%   w is a long column vector.
%   ws is a structure which does not necessarily contain the coefficients
%   OUTPUT
%   ws is the full structure accepted by Waverec4.
%
%   T H   2021
ws.dec = cell(ws.level*15+1,1);
wsz = prod(ws.sizes,2);

ws.dec{1} = reshape(w(1:wsz(1)),ws.sizes(1,:)); % Approximation coefficients
w(1:wsz(1)) = [];

for l = 1:ws.level
    for k = 1:15
        di = 1+15*(l-1)+k; % dec. ind.
        wi = (k-1)*wsz(l);
        ws.dec{di} = reshape(w(wi+1:wi+wsz(l)),ws.sizes(l,:));
    end % k
    w(1:15*wsz(l)) = [];
end % l
end