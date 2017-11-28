SIGNAL =  'Timit' %'Timit' % 'Timit' 'Car'
NOISE =   'Dir_Crown' %'Dir_Crown' %  'Dir_White_Gauss' 'Car'
M = 2;  % No. of Microphones
SNR =    3*ones(M,1); 
snr = 100;

[z,T_sg,vst,T,Tblk,s,v] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g,FLE_l,FLE_r,u,yC,yAS] = algorithm(z,T_sg,vst);

[SNRin,SNRknw,SNRnst,h_known] = blocking(s,v,h,FLE_l,FLE_r,T,Tblk);
