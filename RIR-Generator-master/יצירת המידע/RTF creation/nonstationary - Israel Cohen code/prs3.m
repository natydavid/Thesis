function prs3
% New approach of the system identification (off-line) for stationary noise

[z,Fs,Nbits]=wavread('z');
% x=wavread('x'); ax=wavread('ax'); z=[x ax];
d=wavread('d'); bd=wavread('bd');

leng =20;
%Nfft = 2^nextpow2(2*leng+1);  % 4*leng;
Nfft=2^7;

dM=Nfft/4;
wintype='Hamming';

X=stft(z(:,1),Nfft,dM,1,wintype);
Y=stft(z(:,2),Nfft,dM,1,wintype);
Sxx=abs(X).^2;
Syx=Y.*conj(X);
clear X Y
% smooth spectrograms:
alpha=0.85;
Sxx=filter(1-alpha,[1 -alpha],Sxx,[],2);
Syx=filter(1-alpha,[1 -alpha],Syx,[],2);

D=stft(d,Nfft,dM,1,wintype);
BD=stft(bd,Nfft,dM,1,wintype);
% Suu=mean(abs(D).^2,2);
% Swu=mean(BD.*conj(D),2);
alpha=0.95;
Suu=filter(1-alpha,[1 -alpha],abs(D).^2,[],2);
Swu=filter(1-alpha,[1 -alpha],BD.*conj(D),[],2);
clear D BD

A = ones(Nfft,1);
for w = 1:Nfft/2+1
%     p=find(Sxx(w,:)>2*Suu(w));
    p=find(Sxx(w,:)>2*Suu(w,:));
    if length(p), A(w) =sum((Syx(w,p)-Swu(w,p)).*(Sxx(w,p)-Suu(w,p)))/sum((Sxx(w,p)-Suu(w,p)).^2); end
end
A(Nfft/2+2:Nfft)=conj(A(Nfft/2:-1:2));
a3 = real(ifft(A));
%a3=a3(1:leng);

% compute signal leakage
[x,Fs,Nbits]=wavread('x');
ax=wavread('ax');
S=stft(x,Nfft,dM,1,wintype);
R=stft(ax,Nfft,dM,1,wintype);
L=size(S,2);
R=R-A(1:Nfft/2+1,ones(L,1)).*S;
r3=istft(R,Nfft,dM,1,wintype);

% N=length(x);
% a3x=conv(a3,x); a3x=a3x(1:N);
% r3=ax-a3x;
wavwrite(r3,Fs,Nbits,'r3');

