function prs4
% New approach of the system identification (on-line) for stationary noise

[z,Fs,Nbits]=wavread('z');
% x=wavread('x'); ax=wavread('ax'); z=[x ax];
d=wavread('d'); bd=wavread('bd');

leng =20;
Nfft = 2^nextpow2(2*leng+1);  % 4*leng;
%Nfft=2^7;

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

% compute signal leakage
[x,Fs,Nbits]=wavread('x');
ax=wavread('ax');
S=stft(x,Nfft,dM,1,wintype);
R=stft(ax,Nfft,dM,1,wintype);

mu=0.01;
%mu=0.05;
A = ones(Nfft/2+1,1);
% for l=1:size(S,2)
%     p=find(Sxx(:,l)>2*Suu);
%     if length(p)
%         Err=Syx(p,l)-Swu(p)-A(p).*(Sxx(p,l)-Suu(p));
%         A(p) =A(p)+mu*(Sxx(p,l)-Suu(p)).*conj(Err); 
%     end
%     R(:,l)=R(:,l)-A.*S(:,l);
% end
% r4=istft(R,Nfft,dM,1,wintype);
% wavwrite(r4,Fs,Nbits,'r4');

% mu=0.05;
% for l=1:size(S,2)
%     p=find(Sxx(:,l)>2*Suu);
%     if length(p)
%         Err=Syx(p,l)-Swu(p)-A(p).*(Sxx(p,l)-Suu(p));
%         A(p) =A(p)+mu*(Sxx(p,l)-Suu(p)).*(sign(real(Err))+j*sign(imag(Err))); 
%     end
%     R(:,l)=R(:,l)-A.*S(:,l);
% end
% r4=istft(R,Nfft,dM,1,wintype);
% wavwrite(r4,Fs,Nbits,'r4');

% mu=0.1;
% %s = repmat(initlms(1,mu),Nfft/2+1,1);
% s = repmat(initnlms(1,mu),Nfft/2+1,1);
% % s = repmat(initkalman(1,0.5,2,0.1),Nfft/2+1,1);
% % s = repmat(initrls(1,1,0.9),Nfft/2+1,1);
% for l=1:size(S,2)
%     p=find(Sxx(:,l)>2*Suu);
%     for k=1:length(p)
% %        [y,e,s(p(k))] = adaptlms(Sxx(p(k),l)-Suu(p(k)),Syx(p(k),l)-Swu(p(k)),s(p(k)));
%        [y,e,s(p(k))] = adaptnlms(Sxx(p(k),l)-Suu(p(k)),Syx(p(k),l)-Swu(p(k)),s(p(k)));
% %         [y,e,s(p(k))] = adaptkalman(Sxx(p(k),l)-Suu(p(k)),Syx(p(k),l)-Swu(p(k)),s(p(k)));
% %        [y,e,s(p(k))] = adaptrls(Sxx(p(k),l)-Suu(p(k)),Syx(p(k),l)-Swu(p(k)),s(p(k)));
%     end
%     R(:,l)=R(:,l)-[s.coeffs].' .* S(:,l);
% end
% r4=istft(R,Nfft,dM,1,wintype);
% wavwrite(r4,Fs,Nbits,'r4');

% mu=0.01;
% for l=1:size(S,2)
%     p=find(Sxx(:,l)>2*Suu);
%     if length(p)
%         Err=Syx(p,l)-Swu(p)-A(p).*(Sxx(p,l)-Suu(p));
%         Ex=(Sxx(p,l)-Suu(p)).^2;
%         A(p) =A(p)+mu./Ex.*(Sxx(p,l)-Suu(p)).*Err; 
%     end
%     R(:,l)=R(:,l)-A.*S(:,l);
% end
% r4=istft(R,Nfft,dM,1,wintype);
% wavwrite(r4,Fs,Nbits,'r4');

mu=0.05;
Ex=ones(Nfft/2+1,1);
for l=1:size(S,2)
    p=find(Sxx(:,l)>2*Suu(:,l));
    if length(p)
        Err=Syx(p,l)-Swu(p,l)-A(p).*(Sxx(p,l)-Suu(p,l));
%         Ex(p)=0.7*Ex(p)+0.3*(Sxx(p,l)-Suu(p)).^2;
        Ex(p)=(Sxx(p,l)-Suu(p,l)).^2;
        A(p) =A(p)+mu./Ex(p).*(Sxx(p,l)-Suu(p,l)).*Err; 
    end
    R(:,l)=R(:,l)-A.*S(:,l);
end
r4=istft(R,Nfft,dM,1,wintype);
wavwrite(r4,Fs,Nbits,'r4');



% mu=0.01;
% Yr=Sxx(:,1)-Suu;
% YC_buf=Syx(:,1)-Swu;
% for l=1:size(S,2)
%     Yr_old=Yr;
%     YC_buf_old = YC_buf;
%     Y_old = YC_buf_old - Yr_old.*A;
%     
%     Yr = Sxx(:,l)-Suu; 
%     YC_buf = Syx(:,l)-Swu; 
%     YAS_buf = Yr.*A;
%     Y = YC_buf - YAS_buf;
%     YC_a2=abs(YC_buf).^2;
%     Ya2=abs(Y).^2;
%     Ya2r=abs(Yr).^2;
%     
%     if l==1
%         Ya_avr=abs(Y);
%         beam_avr=YC_a2;
%         ref_avr=beam_avr;
%     else
%         Ya_avr=0.95*Ya_avr+0.05*abs(Y);
%     end
% %     idx=find(abs(S(:,l)).^2<Suu/5);
%     idx=find(Sxx(:,l)<2*Suu);
%     if length(idx)
%         beam_avr(idx)=0.9*beam_avr(idx)+0.1*YC_a2(idx);
%         ref_avr(idx)=0.9*ref_avr(idx)+0.1*Ya2r(idx);
%     end
%     p=find(Sxx(:,l)>2*Suu);
%     if length(p)
%         rel=beam_avr(p)./ref_avr(p);
%         W1_grad=0.5*rel.*mu.*(conj(Yr_old(p)).*( sign(real(Y_old(p)))+j*sign(imag(Y_old(p))))+conj(Yr(p)).*(sign(real(Y(p)))+j*sign(imag(Y(p)) ) ) );
%         A(p) =A(p)+W1_grad./Ya_avr(p); 
%     end
%     R(:,l)=R(:,l)-A.*S(:,l);
% end
% r4=istft(R,Nfft,dM,1,wintype);
% wavwrite(r4,Fs,Nbits,'r4');
