disp('RESULTS_fr_cr_0.mat')

SIGNAL =  'Timit' % 'Timit_br' 
NOISE =  'Dir_Crown'  %  'Dir_Crown' 'Dir_White_Gauss' 'Diffuse' 'Ucorr'
M = 5;  % No. of Microphones
SNR =    0*ones(M,1); 
snr = 100;


[z,T_sg,vst,s] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
ss = s(1,:);
yMX = ncanSH1(yO');
%save RESULTS_fr_cr_-10.mat zz h g yO yMX

% --------------------------------------------------
