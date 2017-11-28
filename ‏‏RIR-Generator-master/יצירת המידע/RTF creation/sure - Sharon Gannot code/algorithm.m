function [yO,h,g,FLE_l,FLE_r,u,yC,yAS] = algorithm(z,T_sg,vst);

[M,N] = size(z);

% ----------------------------
% Algorithm's Parameters
% ----------------------------

mue = 0.05;
eta = 0.05;
lambda = 0.9;
gamma = 1;
FLE_l= 125; FLE_r = 125; % 90
K_l = 125; K_r = 125;


L = 256; % segment length = 2*L
NFFT = 2*L;
JF = NFFT - (FLE_l+FLE_r-1); % Jump

h = zeros(M,FLE_l+FLE_r);
h(1,FLE_l+1) = 1;


g = zeros(M-1,NFFT);
G = fft(g,NFFT,2);



yO = zeros(1,N);
yC = zeros(1,N);
yAS = zeros(1,N);
u = zeros(M-1,N);
yA = zeros(M-1,N);
z_MFBF = zeros(M,N);


%-------------------
% ATF Estimation
%-------------------

% For a while do it outside. Next put into loop !!!!

for m = 2:M
  [H1(m,:),h(m,:)] = est2frq1(z(1,:),z(m,:),T_sg,FLE_l,FLE_r);
end;
H = fft([h(:,FLE_l+1:FLE_l+FLE_r) , zeros(M,NFFT-FLE_l-FLE_r) , h(:,1:FLE_l)],NFFT,2);

Hnorm = sum(abs(H).^2,1);

% -----------------------
% ALGORITHM
% -----------------------

for n = ceil( (NFFT-1)/JF ):fix(N/JF)
  
  T = n*JF-NFFT+1:n*JF; Tc = T(FLE_r:end-FLE_l);

  z_buf = z(:,T);
  Z_buf = fft(z_buf,NFFT,2);
  

  U_buf = Z_buf(2:M,:) - H(2:M,:) .* repmat(Z_buf(1,:),M-1,1); 
  
  u_buf = real(ifft(U_buf,NFFT,2));
  u(:,Tc) = u_buf(:,FLE_r:end-FLE_l);

  Z_MFBF_buf = Z_buf .* conj(H); 

  z_MFBF_buf = real(ifft(Z_MFBF_buf,NFFT,2));
  z_MFBF(:,Tc) = z_MFBF_buf(:,FLE_r:end-FLE_l); 

  YC_buf = sum(Z_MFBF_buf,1) ./ Hnorm; 
  
  yC_buf = real(ifft(YC_buf,NFFT));
  yC(Tc) = yC_buf(FLE_r:end-FLE_l); 

end;

JK = NFFT - (K_l+K_r-1); % Jump
P = zeros(M-1,NFFT);

vss = zeros(ceil(N/JK));

for ii = 1:size(vst,1)
  vss(ceil(vst(ii,1)/JK):ceil(vst(ii,2)/JK)) = 1;
end  


for n = ceil( (NFFT-1)/JK ):fix(N/JK)
  
  T = n*JK-NFFT+1:n*JK;  Tc = T(K_r:end-K_l);

  yC_buf = yC(T);
  U_buf = fft(u(:,T),NFFT,2);

  
  YA_buf = U_buf .* G; 

  yA_buf = real(ifft(YA_buf,NFFT,2));
  yA(:,Tc) = yA_buf(:,K_r:end-K_l); 

  YAS_buf = sum(YA_buf,1);

  yAS_buf = real(ifft(YAS_buf,NFFT));
  yAS(Tc) = yAS_buf(K_r:end-K_l); 
  
  yO_buf = yC_buf - yAS_buf;

  yO(Tc) = yO_buf(K_r:end-K_l);
  YO_buf = fft([zeros(1,K_r-1) , yO_buf(K_r:end-K_l) , zeros(1,K_l)],NFFT);
  
  P = lambda*P + (1-lambda)*abs(U_buf.^2) + eps; 

  TEMP = repmat(YO_buf,M-1,1) .* conj(U_buf) ./ P ;
  Temp = real(ifft(TEMP,NFFT,2));
  TEMP = fft([Temp(:,1:K_r) , zeros(M-1,NFFT-K_l-K_r) , Temp(:,end-K_l+1:end)],NFFT,2);

  G = gamma*G + mue*vss(n)*TEMP;
  
  if ~rem(n,10), 
    disp(['Frame No. ',num2str(n),' Out of ',num2str(fix((N)/JK)),' Frames']);
  end;
  
end

gg = real(ifft(G,NFFT,2));
g = [gg(:,end-K_l+1:end) , gg(:,1:K_r)];


% --------------------------
% End of Algorithm
% --------------------------


