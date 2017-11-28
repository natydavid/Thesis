disp('RESULTS_fr_cr_-10.mat')

SIGNAL =  'Timit' % 'Timit' 'Timit_br'
NOISE =   'Dir_Crown' % 'Dir_White_Gauss''Dir_Crown''Diffuse''Ucorr'
snr = 100;
ts2=8.8e4:9.8e4;
tn=7e4:8.4e4;

for M  = 2:4  % No. of Microphones

  SNR =    -5*ones(M,1); 



  [z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
  [yO,h,g] = algorithm(z,T_sg,vst);
  zz = z(1,:);
  yMX = ncanSH1(yO');
  zMX = ncanSH1(zz');
  
  SNR_z(M) = nr(zz,ts2,tn)
  SNR_o(M) = nr(yO,ts2,tn)
  SNR_mx(M) = nr(yMX,ts2,tn)
  
end;

%save RESULTS_fr_cr_-10.mat zz h g yO yMX

% --------------------------------------------------
