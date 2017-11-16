
disp('RESULTS_fr_cr_-10.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'Crown' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    -10*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_cr_-10.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_cr_-7.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'Crown' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    -7*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_cr_-7.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_cr_-5.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'Crown' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    -5*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_cr_-5.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_cr_-3.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'Crown' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    -3*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_cr_-3.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_cr_0.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'Crown' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    0*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_cr_0.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_cr_3.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'Crown' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    3*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_cr_3.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_cr_6.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'Dir_Crown' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    5*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_cr_6.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_cr_10.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'Crown' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    10*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_cr_10.mat zz h g yO yMX

% -------------------------------------------
% -------------------------------------------
% -------------------------------------------
% -------------------------------------------

disp('RESULTS_fr_wt_-10.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'White_Gauss' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    -10*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_wt_-10.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_wt_-7.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'White_Gauss' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    -7*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_wt_-7.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_wt_-5.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'White_Gauss' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    -5*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_wt_-5.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_wt_-3.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'White_Gauss' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    -3*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_wt_-3.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_wt_0.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'White_Gauss' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    0*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_wt_0.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_wt_3.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'White_Gauss' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    3*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_wt_3.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_wt_5.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'White_Gauss' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    5*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_wt_5.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_wt_7.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'White_Gauss' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    7*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_wt_7.mat zz h g yO yMX

% --------------------------------------------------
disp('RESULTS_fr_wt_10.mat')

SIGNAL =  'Timit' % 'Timit' 
NOISE =   'White_Gauss' %  'Crown' 'White_Gauss'
M = 5;  % No. of Microphones
SNR =    10*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
yMX = ncanSH1(yO');
save RESULTS_fr_wt_10.mat zz h g yO yMX
