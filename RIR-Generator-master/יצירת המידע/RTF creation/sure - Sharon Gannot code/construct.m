function [Z,T_sg,vst,T,Tblk,S,V] = construct(SIGNAL,NOISE,SNR,snr,M);

% M = length(SNR);

% -----------------------
% DESIRED SIGNAL
% -----------------------

switch SIGNAL
    
case 'Car'
    
    S(1,:) = wavread('x1')';
    S(2,:) = wavread('x2')';
    S(3,:) = wavread('x3')';
    S(4,:) = wavread('x4')';
    
    N = length(S);
    S = S(1:M,1:N);
  
    %     T_sg = [.28e4,.38e4;
    %         1.21e4,1.31e4;
    %         2.94e4,3.04e4;
    %         3.85e4,3.95e4;
    %         4.68e4,4.78e4;
    %         6.32e4,6.42e4;
    %         7.22e4,7.32e4;
    %         7.33e4,7.43e4];
    
    T_sg = [2801,2900;
        2901,3000;
        3001,3100;
        3101,3200;
        3201,3300;
        3301,3400;
        3401,3500;
        3501,3600;
        3601,3700;
        3701,3800];
    %T = 5.55e4:5.91e4;
    
    T = 0.25e4:0.45e4;%7.33e4:7.43e4;%
    Tblk = 7.33e4:7.43e4;%7.22e4,7.32e4 % 0.25e4:0.45e4;
    
    
    vst = [0.70e4 , 1.00e4;
        1.55e4 , 2.00e4;
        2.30e5 , 2.65e4;
        3.25e4 , 3.55e4;
        4.15e4 , 4.40e4;
        5.08e4 , 5.24e4;
        5.90e4 , 6.08e4;
        6.75e4 , 7.15e4;
        7.70e4 , 7.9e4];
    

case 'Timit'
    
    S(1,:) = rd('ti1.2','l')';
    S(2,:) = rd('ti1.3','l')';
    S(3,:) = rd('ti1.4','l')';
    S(4,:) = rd('ti1.5','l')';
    S(5,:) = rd('ti1.6','l')';
    
    
    N = 18e4;
    S = S(1:M,1:N);
    
    T_sg = [5.6e4,5.7e4;
        5.7e4,5.8e4;
        5.8e4,5.9e4;
        8.6e4,8.7e4;
        8.7e4,8.8e4;
        8.8e4,8.9e4;
        8.9e4,9.0e4;
        9.0e4,9.1e4;
        9.1e4,9.2e4;
        9.2e4,9.3e4;
        9.3e4,9.4e4;
        9.4e4,9.5e4;
        9.5e4,9.6e4];
    
    %T = 5.55e4:5.91e4;
    
    T = 8.6e4:9.8e4;
    Tblk = 8.6e4:9.8e4;
    
    
    vst = [1      , 49970;
        69958  , 81793;
        103885 , 115983;
        137812 , 151751];
    % vst(1:49970) = 1; vst(69958:81793) = 1;
    % vst(103885:115983) = 1; vst(137812:151751) = 1; 
    
case 'Timit_br'
    
    S(1,:) = rd('ti2br.2','l')';
    S(2,:) = rd('ti2br.3','l')';
    S(3,:) = rd('ti2br.4','l')';
    S(4,:) = rd('ti2br.5','l')';
    S(5,:) = rd('ti2br.6','l')';
    
    
    N = 15e4;
    S = S(1:M,1:N);
    
    T_sg = [1.9e4,2.0e4;
        2.0e4,2.1e4;
        2.1e4,2.2e4;
        5.3e4,5.4e4;
        5.4e4,5.5e4;
        5.5e4,5.6e4;
        5.7e4,5.8e4;
        5.9e4,6.0e4;
        6.0e4,6.1e4
        6.1e4,6.2e4
        6.2e4,6.3e4
        6.3e4,6.4e4
        6.4e4,6.5e4];
    
    
    
    vst = [1      , 1.0e4;
        3.8e4  , 4.8e4;
        7.0e4  , 8.2e4;
        10.6e4 , 11.6e4];
    
otherwise
    
    error('No other signal available')
    
end;

% -----------------------
% INTERFERENCE SIGNAL
% -----------------------

switch NOISE
    
case 'Car'
    
    VV(1,:) = wavread('n1')';
    VV(2,:) = wavread('n2')';
    VV(3,:) = wavread('n3')';
    VV(4,:) = wavread('n4')';
    
    VV = VV(1:M,1:N);
    
case 'Dir_White_Gauss'
    
    VV(1,:) = rd('wt.2','l')';
    VV(2,:) = rd('wt.3','l')';
    VV(3,:) = rd('wt.4','l')';
    VV(4,:) = rd('wt.5','l')';
    VV(5,:) = rd('wt.6','l')';
    VV(6,:) = rd('wt.7','l')';
    
    n = 1.1e5;
    VV = VV(1:M,n+1:n+N);
    
case 'Dir_Crown'
    
    VV(1,:) = rd('cr.2','l')';
    VV(2,:) = rd('cr.3','l')';
    VV(3,:) = rd('cr.4','l')';
    VV(4,:) = rd('cr.5','l')';
    VV(5,:) = rd('cr.6','l')';
    VV(6,:) = rd('cr.7','l')';
    
    n = 1.1e5;
    VV = VV(1:M,n+1:n+N);
    
case 'Diffuse'
    
    VV(1,:) = rd('df.2','l')';
    VV(2,:) = rd('df.3','l')';
    VV(3,:) = rd('df.4','l')';
    VV(4,:) = rd('df.5','l')';
    VV(5,:) = rd('df.6','l')';
    
    n = 1e5;
    VV = VV(1:M,n+1:n+N);
    
    [B,A] = butter(5,.2,'high');
    
    for m = 1:M
        VV(m,:) = filter(B,A,VV(m,:));
    end;
    
case 'Ucorr'
    
    randn('state',0)
    VV = 1000*randn(M,N);
    
    [B,A] = butter(5,.2,'high');
    
    for m = 1:M
        VV(m,:) = filter(B,A,VV(m,:));
    end;
    
otherwise
    
    error('No other noise source available')
    
end;


% -----------------------
% Signals' construction
% -----------------------

gv = sqrt(var(S(1,:))*10^(-snr/10));
NO = gv*randn(size(S));


for m = 1:M
    
    GV(m) = sqrt(var(S(m,:))/var(VV(m,:))*10^(-SNR(m)/10));
    
end;


V = (GV'*ones(1,length(S))).*VV;
Z = S + V + NO;



% ----- End of construct ---------

