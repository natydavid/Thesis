function [Z,T_sg,T,Tblk,vst,S,V] = construct(SIGNAL,NOISE,SNR,snr,M);

% M = length(SNR);

% -----------------------
% DESIRED SIGNAL
% -----------------------

switch SIGNAL
 
case 'Image'
    
    s = [zeros(1e4,1) ; rd('eyes.8k','b') ; zeros(1e4,1) ; rd('doors.8k','b') ; zeros(1e4,1)];
    N = length(s); s = s(1:N);

    T_sg =...
        [1.2e4,1.3e4;
         1.3e4,1.4e4;
         1.8e4,1.9e4;
         1.9e4,2.0e4;
         2.3e4,2.4e4;
         2.4e4,2.5e4;
         3.8e4,3.9e4;
         3.9e4,4.0e4;
         4.2e4,4.3e4];
    
    T = 1.2e4:2.4e4;
    
    Tblk = 2.7e4:3.5e4;
    
    
    vst = [1   , .9e4;
           2.7e4  , 3.5e4;
           5.3e4  , 6.1e4];
    
    % Artificial room:

    dim = [1 1 1];
    src = [.1 .9 .1];
    Tr = 0.1;
    L = 1024;
    
    h = simulate([.1 .1 .1],src,dim,Tr,L); a(1,:) = h.h';
    h = simulate([.3 .1 .1],src,dim,Tr,L); a(2,:) = h.h';
    h = simulate([.5 .1 .1],src,dim,Tr,L); a(3,:) = h.h';
    h = simulate([.7 .1 .1],src,dim,Tr,L); a(4,:) = h.h';
    h = simulate([.9 .1 .1],src,dim,Tr,L); a(5,:) = h.h';
    
    for m = 1:M    
        S(m,:) = filter(a(m,:),1,s');
    end;
    
    
    
    
case 'Speech'
    
    s = [zeros(1e4,1) ; rd('Agnes.8k','b')];
    N = 60000; s = s(1:N);
    T_sg = round(...
        [.54 .59
        .63 .68
        2.45 2.5
        2.55 2.6
        2.65 2.7
        2.75 2.8
        2.81 2.86
        3.2  3.25]*1e4)+1e4;
    %5.25 5.3
    %5.5  5.55]*1e4);
    T = 2.5e4:3.5e4;
    Tblk = 2.5e4:3.5e4;%15001:30000;
    
    
    a(1,:) = [1 0 0 0 0 0 0 0 0 0];
    a(2,:) = [0 0 1 0 0 0 0 0 0 0];
    a(3,:) = [0 0 0 0 1 0 0 0 0 0];
    a(4,:) = [0 0 0 0 0 0 1 0 0 0];
    a(5,:) = [0 0 0 0 0 0 0 0 1 0];
    
    for m = 1:M
        
        S(m,:) = filter(a(m,:),1,s');
        
    end;
    
    vst = [1      , 49970];
    
    
case 'Car'
    
    S(1,:) = wavread('x1')';
    S(2,:) = wavread('x2')';
    S(3,:) = wavread('x3')';
    S(4,:) = wavread('x4')';
    
    
    N = length(S);
    S = S(1:M,1:N);
    
    
%     T_sg = [2801,2900;
%         2901,3000;
%         3001,3100;
%         3101,3200;
%         3201,3300;
%         3301,3400;
%         3401,3500;
%         3501,3600;
%         3601,3700;
%         3701,3800];
% 
  
    T_sg =...
        [2801,3300;
        3301,3800;
        3801,4300;
        12001,12500;
        12501,13000;
        46501, 47000;
        47001, 47500;
        47501, 48000];

  

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
  
    case 'Echo'
    
    S(1,:) = wavread('xE1')';
    S(2,:) = wavread('xE2')';
    S(3,:) = wavread('xE3')';
    S(4,:) = wavread('xE4')';
    
    
    N = length(S);
    S = S(1:M,1:N);
  
    T_sg =...
        [53201  ,  53700;
         53701  ,  54200;
         54201  ,  54700;
         54701  ,  55200;
         135001 , 135500;
         135501 , 136000;
         144001 , 144500;
         144501 , 145000;
         145001 , 145500;
         145501 , 146000
     ];

  

    T = [5.32e4:5.52e4,1.35e5:1.36e5,1.44e5:1.46e5];
    Tblk = 1.05e5:1.07e5;
    
    
    vst = [...
            5.80e4 , 6.30e4;
            7.05e4 , 7.50e4;
            7.90e4 , 8.20e4;
            8.95e4 , 9.20e4;
            0.99e5 , 1.025e5;
            1.105e5 , 1.125e5;
            1.205e5 , 1.225e5;
            1.305e5 , 1.335e5;
            1.402e5 , 1.422e5;...
        ];

 case 'Car_Spk'
    
    ss = wavread('x1')';
    load Left;
    
    S(1,:) = filter(Ll(1,:),1,ss);
    S(2,:) = filter(Ll(2,:),1,ss);
    S(3,:) = filter(Ll(3,:),1,ss);
    S(4,:) = filter(Ll(4,:),1,ss);
    S(5,:) = filter(Ll(5,:),1,ss);
    
    
    N = length(S);
    S = S(1:M,1:N);
    
    T_sg =...
        [2801,3300;
        3301,3800;
        3801,4300;
        12001,12500;
        12501,13000;
        46501, 47000;
        47001, 47500;
        47501, 48000];

%     T_sg = [2801,2900;
%         2901,3000;
%         3001,3100;
%         3101,3200;
%         3201,3300;
%         3301,3400;
%         3401,3500;
%         3501,3600;
%         3601,3700;
%         3701,3800;
%         3801,3900;
%         3901,4000;
%         4001,4100];
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
    Ni = 4e4;
    S = S(1:M,Ni+1:Ni+N);
    
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
    
    T = 5.55e4:5.91e4;
    
    %T = 8.6e4:9.8e4;
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
    
case 'Image'
    
    % Artificial room:

    dim = [1 1 1];
    src = [.9 .9 .1];
    Tr = 0.1;
    L = 1024;

    randn('state',0)
    v = 1000*randn(N,1);
    
    h = simulate([.1 .1 .1],src,dim,Tr,L); b(1,:) = h.h';
    h = simulate([.3 .1 .1],src,dim,Tr,L); b(2,:) = h.h';
    h = simulate([.5 .1 .1],src,dim,Tr,L); b(3,:) = h.h';
    h = simulate([.7 .1 .1],src,dim,Tr,L); b(4,:) = h.h';
    h = simulate([.9 .1 .1],src,dim,Tr,L); b(5,:) = h.h';

    for m = 1:M    
        VV(m,:) = filter(b(m,:),1,v');
    end;

case 'Car'
    
    VV(1,:) = wavread('n1')';
    VV(2,:) = wavread('n2')';
    VV(3,:) = wavread('n3')';
    VV(4,:) = wavread('n4')';
    
    VV = VV(1:M,1:N);
    
case 'Echo'
    
    VV(1,:) = wavread('nE1')';
    VV(2,:) = wavread('nE2')';
    VV(3,:) = wavread('nE3')';
    VV(4,:) = wavread('nE4')';
    
    VV = VV(1:M,1:N);
    
case 'Car_Spk'
    
    vv = wavread('n1')';
    load Right;
    
    VV(1,:) = filter(Rr(1,:),1,vv);
    VV(2,:) = filter(Rr(2,:),1,vv);
    VV(3,:) = filter(Rr(3,:),1,vv);
    VV(4,:) = filter(Rr(4,:),1,vv);
    VV(5,:) = filter(Rr(5,:),1,vv);
        
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
    
case 'Diff_NS'
    
    VV(1,:) = rd('dfns.2','l')';
    VV(2,:) = rd('dfns.3','l')';
    VV(3,:) = rd('dfns.4','l')';
    VV(4,:) = rd('dfns.5','l')';
    VV(5,:) = rd('dfns.6','l')';
    
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

