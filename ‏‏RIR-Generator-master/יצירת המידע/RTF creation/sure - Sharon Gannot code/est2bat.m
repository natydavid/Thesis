function [a] = est2bat(x,y,L1,L2);
%
L22 = L2-1;
L = max(L1,L22);

rxx = xcorr(x,L1+L22,'none');
%ryx = xcorr(x,y,L,'none');
ryx = xcorr(y,x,L,'none');
Rxx = toeplitz(rxx(L1+L22+1:2*(L1+L22)+1));
Ryx = [ryx(L-L1+1:L) ; ryx(L+1:L+L22+1)];
a = (Rxx \ Ryx)';

