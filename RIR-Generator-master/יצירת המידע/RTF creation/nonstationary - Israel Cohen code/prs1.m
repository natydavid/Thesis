function prs1
% create clean and noisy signals

% read clean signal
[x,Fs,Nbits]=wavread('draw_8');
x=[zeros(0.6*Fs,1); 0.25*x];      % start with 0.6 s of silence
N=length(x);

%%%%%%%%%%%%%
%x=0.1*randn(N,1);
%%%%%%%%%%%%%

wavwrite(x,Fs,Nbits,'x');
Ex=mean(x.^2);

% create noise
new_noise=0;
if new_noise    % create a new noise signal
    %  Ed=2*Ex;
    Ed=0.5*Ex;
    %     Ad=[linspace(0.5,2,ceil(N/2)) linspace(2,0.5,floor(N/2))]';
    %     Ad=[linspace(0.5,2.5,ceil(N/4)) linspace(2.5,0.5,floor(N/4))]'; Ad=[0; Ad; Ad];
    %    Ad=3*(1+cos(2*pi*(0:N-1)'/Fs));
    Ad=[linspace(1,4,ceil(N/2)) linspace(4,1,floor(N/2))]';
    d=Ed^0.5*Ad.*randn(N,1);
    wavwrite(d,Fs,Nbits,'d');
else
    d=wavread('d');
end

SegSNR=segsnr(x(:,1),d(:,1))

% create noisy signals
%a=[1 0.1 0.2]';
%a=[0 0 1 0 0 0.7]';
%a=[0 0 0 0 0 0 0.1 0.2 0.3]';
a=[0 0 0 0 0 0 1 -0.5 0.25]';
%a=a/norm(a);
b=[-1 -0.5 0.1]';
%b=b/norm(b);
ax=conv(a,x); ax=ax(1:N);
bd=conv(b,d); bd=bd(1:N);
wavwrite(ax,Fs,Nbits,'ax');
wavwrite(bd,Fs,Nbits,'bd');
z=[x+d ax+bd];
wavwrite(z,Fs,Nbits,'z');
wavwrite(z(:,1),Fs,Nbits,'z1');
wavwrite(z(:,2),Fs,Nbits,'z2');

function d=segsnr(x,d)
% Segmental SNR
M=256;
Mo=0.5*M;
M21=M/2+1;
X=stft(x,M,M-Mo,1);
D=stft(d,M,M-Mo,1);
X=sum(abs(X(3:M21,:)).^2);
D=sum(abs(D(3:M21,:)).^2);
d=repmat(-15,size(X));
%p=find(X>max(X)*1e-6);
p=find(X);
d=10*log10(X(p)./D(p));
d(find(d>35))=35;
d=mean(d);
