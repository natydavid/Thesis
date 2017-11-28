SIGNAL =  'Timit' % 'Timit' 
NOISE =   'Dir_Crown' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    -3*ones(M,1); 
snr = 100;


[z,T_sg,T,Tblk,vst,s,v] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g,FLE_l,FLE_r] = algorithm(z,T_sg,vst);

[SNRin,SNRknw,SNRnst] = blocking(s,v,h,FLE_l,FLE_r,T,Tblk);

tsn = 8.9e4:9.6e4;
tsn = 8.9e4:9.6e4;

nr(z(1,:),tsn,tn)
nr(yO,tsn,tn)

DIST = dist(z(1,:),yO,s(1,:),tsn)
