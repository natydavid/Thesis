function prs2
% system identification according to Shalvi-Weinstein

[z,Fs,Nbits]=wavread('z');
% x=wavread('x'); ax=wavread('ax'); z=[x ax];
x=wavread('x'); d=wavread('d');

% segments containing speech (of length 128 samples)
% T_sg=[17750 19579 20893 22902 23727 25195 28055 31842 33954 39803]'-1.4*Fs;
N=size(z,1);
T_sg=(1:128:N-128)';
T_sg =[T_sg  T_sg+127];

leng =20;
%win = bartlett(2*leng+1);
BL = size(T_sg,1);  % number of segments
Nfft = 2^nextpow2(2*leng+1);  % 4*leng;
Sxx = zeros(BL,Nfft);
Syx = zeros(BL,Nfft);
maxlag = leng;
BLn=0;
for bl = 1 : BL
    T = T_sg(bl,1):T_sg(bl,2);
    SNR=sum(x(T).^2)/sum(d(T).^2);
    if SNR>1
        BLn=BLn+1;
        
        x1 = z(T,1);
        y1 = z(T,2);
        
        %     r_xx = xcorr(x1,maxlag,'unbiased').*win;
        %     r_yx = xcorr(y1,x1,maxlag,'unbiased').*win;
        %     r_xx = xcorr(x1,maxlag,'none').*win;
        %     r_yx = xcorr(y1,x1,maxlag,'none').*win;
        r_xx = xcorr(x1,maxlag,'unbiased');
        r_yx = xcorr(y1,x1,maxlag,'unbiased');
        K_xx = zeros(1,Nfft); 
        K_xx(1:maxlag+1) = r_xx(maxlag+1:2*maxlag+1); 
        K_xx(Nfft-maxlag+1:Nfft) = r_xx(1:maxlag);
        K_yx = zeros(1,Nfft); 
        K_yx(1:maxlag+1) = r_yx(maxlag+1:2*maxlag+1); 
        K_yx(Nfft-maxlag+1:Nfft) = r_yx(1:maxlag);
        
        Sxx(BLn,:) = abs(fft(K_xx,Nfft));
        Syx(BLn,:) = fft(K_yx,Nfft);
    end
end;
Sxx=Sxx(1:BLn,:);
Syx=Syx(1:BLn,:);

% weighted algorithm
A = zeros(Nfft,1);
for w = 1:Nfft 
    A(w)=(sum(Syx(:,w))*sum(1./Sxx(:,w))-sum(Syx(:,w)./Sxx(:,w)))/(sum(Sxx(:,w))*sum(1./Sxx(:,w))-1);
end
a2 = real(ifft(A));

% nonweighted algorithm
A = zeros(Nfft,1);
for w = 1:Nfft 
    A(w)=(sum(Syx(:,w).*Sxx(:,w))-sum(Syx(:,w))*sum(Sxx(:,w)))/(sum(Sxx(:,w).^2)-sum(Sxx(:,w))^2);
end
a1 = real(ifft(A));

% compute signal leakage
[x,Fs,Nbits]=wavread('x');
ax=wavread('ax');

N=length(x);
a1x=conv(a1,x); a1x=a1x(1:N);
r1=ax-a1x;
wavwrite(r1,Fs,Nbits,'r2a');

a2x=conv(a2,x); a2x=a2x(1:N);
r2=ax-a2x;
wavwrite(r2,Fs,Nbits,'r2b');