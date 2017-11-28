function [H,h,Sxx,Syx,Svx] = est2frq1(x,y,T_sg,FLE_l,FLE_r);
%
% LS fitting for the system using non-statinarity (using Matlab LS)

%
% Estimates auto and cross-spectra of x and y within sub-periods specified
% by  'periods'.
% Usage: [Sxx,Syx]=CalcSpec(x,y,periods,leng);
%
% Input: x, y - signals
%        periods - length of each segment (the spectrum is calculated
%        seperatly in each segment
%        leng -  2*leng+1 is the window length in Blackman- Tukey method  for spectral-estimation
%        (the window type is specified inside this program - Hamming).
%

leng = 1*(FLE_l + FLE_r);

win = hamming(2*leng+1)';

BL = length(T_sg);


Nfft = 2^nextpow2(2*leng+1); % 4*leng;

Sxx = zeros(BL,Nfft);
Syx = zeros(BL,Nfft);

for bl = 1 : BL


  T = T_sg(bl,1):T_sg(bl,2);
  

  x1 = x(T);
  y1 = y(T);

  
  % Blackman-Tukey


  maxlag = leng;

  r_xx = xcorr(x(T),maxlag,'none').*win;
  r_yx = xcorr(y(T),x(T),maxlag,'none').*win;%xcorr(x(T),y(T),maxlag,'none').*win;
  K_xx = zeros(1,Nfft); 
  K_xx(1:maxlag+1) = r_xx(maxlag+1:2*maxlag+1); 
  K_xx(Nfft-maxlag+1:Nfft) = r_xx(1:maxlag);
  K_yx = zeros(1,Nfft); 
  K_yx(1:maxlag+1) = r_yx(maxlag+1:2*maxlag+1); 
  K_yx(Nfft-maxlag+1:Nfft) = r_yx(1:maxlag);
 
  Sxx(bl,:) = fft(K_xx,Nfft);
  Syx(bl,:) = fft(K_yx,Nfft);


end;

teta = zeros(2,Nfft/2+1);

for w = 19:200 %225 Nfft/2+1 % 6:120% 16:64 % looping for each frequancy

  [teta(:,w)] = [Sxx(:,w) , ones(size(Sxx(:,w)))] \ Syx(:,w);

end

H   = [teta(1,1:Nfft/2+1) conj(teta(1,Nfft/2:-1:2))]; % teta(1,:);%
Svx = [teta(2,1:Nfft/2+1) conj(teta(2,Nfft/2:-1:2))]; % teta(2,:);%

hh = ifft(H,[],2);


h = real([hh(Nfft-FLE_l+1:Nfft) hh(1:FLE_r)]); % synchronized!