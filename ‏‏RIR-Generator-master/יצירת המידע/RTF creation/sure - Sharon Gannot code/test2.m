%disp('RESULTS_fr_cr_-10.mat')

SIGNAL =  'Timit' % 'Timit' 'Timit_br'
NOISE =   'Ucorr' % 'Dir_White_Gauss''Dir_Crown' 'Diffuse''Ucorr'
snr = 100;
ts2=8.8e4:9.8e4;
tn=7e4:8.4e4;

M = 5;

SN = [-10 -7 -5 -3 0 3 6 10];

for m = 1:8 
  
  SN(m)
  SNR = SN(m)*ones(M,1); 



  [z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
  [yO,h,g] = algorithm(z,T_sg,vst);
  zz = z(1,:);
  yMX = ncanSH1(yO');
  zMX = ncanSH1(zz');
  
  SNR_z(m) = nr(zz,ts2,tn)
  SNR_o(m) = nr(yO,ts2,tn)
  SNR_mx(m) = nr(yMX,ts2,tn)
  
end;



%save RESULTS_fr_cr_-10.mat zz h g yO yMX

% --------------------------------------------------
