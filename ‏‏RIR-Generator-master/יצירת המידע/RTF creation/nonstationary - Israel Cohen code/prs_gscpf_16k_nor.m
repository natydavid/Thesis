function prs_gscpf_16k_nor
% normal form of prs_gscpf_16k

parms={[8e3 256 192] ...
        [1 0.9 8 15 1.67 1.66 4.6 3 0.85] ...
        [0.7 1 15 10e3 500 0.005 -5 -10 -5 -10 -5 -10 10 0 0.998] ...
        [0.96 -20] ...
        0};

dir_noise='D:\Data\DA310\';
dir_proc='D:\Data\DA310_processed\';

% signal_name='tmp';
% fin={[dir_noise signal_name '_13']};
signal_name='usb3';
fin={[dir_noise signal_name]};

fbeam=[dir_proc signal_name '_beam'];
fref=[dir_proc signal_name '_ref'];
fout=[dir_proc signal_name '_gscpf_nor'];
gscpf(fin,fbeam,fout,1,parms)

return

parms={[16e3 512 384] ...
        [1 0.9 8 15 1.67 1.66 4.6 3 0.85] ...
        [0.7 1 15 10e3 500 0.005 -5 -10 -5 -10 -5 -10 10 0 0.998] ...
        [0.96 -20] ...
        0};

dir_noise='D:\Data\DA310_2ch\';
dir_proc='D:\Data\DA310_2ch_processed_nor\';

signal_name='Training';
fin={[dir_noise signal_name '_13']};
fbeam=[dir_proc signal_name '_beam'];
fout=[dir_proc signal_name '_gscpf'];
gscpf(fin,fbeam,fout,1,parms)

signal_name='Philos';
fin={[dir_noise signal_name '_13_n1']};
fbeam=[dir_proc signal_name '_n1_beam'];
fout=[dir_proc signal_name '_n1_gscpf'];
gscpf(fin,fbeam,fout,1,parms)

signal_name='Philos';
fin={[dir_noise signal_name '_13_n2']};
fbeam=[dir_proc signal_name '_n2_beam'];
fout=[dir_proc signal_name '_n2_gscpf'];
gscpf(fin,fbeam,fout,1,parms)


function gscpf(fin,fbeam,fout,noise_factor,parms,fref,out_st)
% gscpf : GSC beamforming with Post-Filtering 
% ***************************************************************@
% Inputs: 
%    fin,  cell of strings: input file names (fin{1}.wav; fin{2}.wav;... )
%    fout, output file name (fout.wav)
%    noise_factor, noise overestimation factor
%    parms, parameters
%    fref, output reference file name (fref.wav)
%    out_st, flag whether to produce mono output (out_st=0) or stereo output with left=right  (out_st=1)
% Output:
%    out, samples of the output file (post filtered beam)
% Usage: 
%    out=gscpf(fin,fout,noise_factor,parms);
%    gscpf(fin);
% Defaults:
%    fout= 'y_gscpf';
%    noise_factor=1;
%    out_st=0;

% Copyright (c) 2002. Dr Israel Cohen. Lamar Signal Processing.
% All rights reserved. Created  10/9/02.
% ***************************************************************@

disp('Microphone Array Beamforming and Post-Filtering.')
disp('Copyright (c) 2002. Dr Israel Cohen.')

nin=nargin;

if nin<7
    out_st=0';
end
if nin<3
    fout='y_gscpf';
end
if nin<4
    noise_factor=1;
end
if nin<5		% set default parameters
    % 1) Parameters of Short Time Fourier Analysis:
    Fs_ref=8e3;		% 1.1) Reference Sampling frequency
    M_ref=256;		% 1.2) Size of analysis window
    Mo_ref=0.75*M_ref;	% 1.3) Number of overlapping samples in consecutive frames
    
    % 2) Parameters of Noise Spectrum Estimate
    w=2;			% 2.1)  Size of frequency smoothing window function=2*w+1
    alpha_s_ref=0.8;	% 2.2)  Recursive averaging parameter for the smoothing operation
    Nwin=8; 	% 2.3)  Resolution of local minima search
    Vwin=15;	
    delta_s=1.83;		
    Bmin=1.98;
    delta_y=4.6;		% 2.4)  Local minimum factor
    delta_yt=3;	
    alpha_d_ref=0.85;	% 2.7)  Recursive averaging parameter for the noise
    
    % 3) Parameters of a Priori Probability for Signal-Absence Estimate
    alpha_xi_ref=0.7;	% 3.1) Recursive averaging parameter
    w_xi_local=1; 	% 3.2) Size of frequency local smoothing window function 
    w_xi_global=15; 	% 3.3) Size of frequency local smoothing window function 
    f_u=10e3; 		% 3.4) Upper frequency threshold for global decision
    f_l=500; 		% 3.5) Lower frequency threshold for global decision
    P_min=0.05; 		% 3.6) Lower bound constraint
    xi_lu_dB=-5; 	% 3.7) Upper threshold for local decision
    xi_ll_dB=-10; 	% 3.8) Lower threshold for local decision
    xi_gu_dB=-5; 	% 3.9) Upper threshold for global decision
    xi_gl_dB=-10; 	% 3.10) Lower threshold for global decision
    xi_fu_dB=-5; 	% 3.11) Upper threshold for local decision
    xi_fl_dB=-10; 	% 3.12) Lower threshold for local decision
    xi_mu_dB=10; 	% 3.13) Upper threshold for xi_m
    xi_ml_dB=0; 		% 3.14) Lower threshold for xi_m
    q_max=0.95; 		% 3.15) Upper limit constraint
    
    % 4) Parameters of "Decision-Directed" a Priori SNR Estimate
    alpha_eta_ref=0.92;	% 4.1) Recursive averaging parameter
    eta_min_dB=-20;	% 4.2) Lower limit constraint
else
    % 1) Parameters of Short Time Fourier Analysis:
    Fs_ref=parms{1}(1);	% 1.1) Reference Sampling frequency
    M_ref=parms{1}(2);	% 1.2) Size of analysis window
    Mo_ref=parms{1}(3);	% 1.3) Number of overlapping samples in consecutive frames
    
    % 2) Parameters of Noise Spectrum Estimate
    w=parms{2}(1);			% 2.1)  Size of frequency smoothing window function=2*w+1
    alpha_s_ref=parms{2}(2);	% 2.2)  Recursive averaging parameter for the smoothing operation
    Nwin=parms{2}(3); 	% 2.3)  Resolution of local minima search
    Vwin=parms{2}(4);	
    delta_s=parms{2}(5);		% 2.4)  Local minimum factor
    Bmin=parms{2}(6);
    delta_y=parms{2}(7);		% 2.4)  Local minimum factor
    delta_yt=parms{2}(8);	
    alpha_d_ref=parms{2}(9);	% 2.7)  Recursive averaging parameter for the noise
    
    % 3) Parameters of a Priori Probability for Signal-Absence Estimate
    alpha_xi_ref=parms{3}(1);	% 3.1) Recursive averaging parameter
    w_xi_local=parms{3}(2); 	% 3.2) Size of frequency local smoothing window function 
    w_xi_global=parms{3}(3); 	% 3.3) Size of frequency local smoothing window function 
    f_u=parms{3}(4); 		% 3.4) Upper frequency threshold for global decision
    f_l=parms{3}(5); 		% 3.5) Lower frequency threshold for global decision
    P_min=parms{3}(6); 		% 3.6) Lower bound constraint
    xi_lu_dB=parms{3}(7); 	% 3.7) Upper threshold for local decision
    xi_ll_dB=parms{3}(8); 	% 3.8) Lower threshold for local decision
    xi_gu_dB=parms{3}(9); 	% 3.9) Upper threshold for global decision
    xi_gl_dB=parms{3}(10); 	% 3.10) Lower threshold for global decision
    xi_fu_dB=parms{3}(11); 	% 3.11) Upper threshold for local decision
    xi_fl_dB=parms{3}(12); 	% 3.12) Lower threshold for local decision
    xi_mu_dB=parms{3}(13); 	% 3.13) Upper threshold for xi_m
    xi_ml_dB=parms{3}(14);	% 3.14) Lower threshold for xi_m
    q_max=parms{3}(15); 		% 3.15) Upper limit constraint
    
    % 4) Parameters of "Decision-Directed" a Priori SNR Estimate
    alpha_eta_ref=parms{4}(1);	% 4.1) Recursive averaging parameter
    eta_min_dB=parms{4}(2);	% 4.2) Lower limit constraint
end
Bbeta=1.4685;

% Read input data 
[N,Fs,NBITS]=wavread(fin{1},0);    % read the size of the data, Fs and NBITS
Nfin=length(fin); % Number of input files
Nch_fin=N(2);   % Number of channels per input file
Nch=Nch_fin*Nfin;    % Number of channels
N=N(1);

% Adjust parameters according to the actual sampling frequency
if Fs~=Fs_ref
    M=2^round(log2(Fs/Fs_ref*M_ref));
    Mo=Mo_ref/M_ref*M;
    alpha_s=alpha_s_ref^(M_ref/M*Fs/Fs_ref);
    alpha_d=alpha_d_ref^(M_ref/M*Fs/Fs_ref);
    alpha_eta=alpha_eta_ref^(M_ref/M*Fs/Fs_ref);
    alpha_xi=alpha_xi_ref^(M_ref/M*Fs/Fs_ref);
else
    M=M_ref;
    Mo=Mo_ref;
    alpha_s=alpha_s_ref;
    alpha_d=alpha_d_ref;
    alpha_eta=alpha_eta_ref;
    alpha_xi=alpha_xi_ref;
end
eta_min=10^(eta_min_dB/10);
G_f=eta_min^0.5;	   % Gain floor
Mno=M-Mo;
Nframes=fix((N-Mo)/(M-Mo));   %  number of frames 
M21=M/2+1;

% window function
win=hamming(M);
% find a normalization factor for the window
win2=win.^2;
W0=win2(1:Mno);
for k=Mno:Mno:M-1
    swin2=lnshift(win2,k);
    W0=W0+swin2(1:Mno);
end
W0=mean(W0)^0.5;
win=win/W0;
Cwin=sum(win.^2)^0.5;
win=win/Cwin;

out=zeros(M,1); % % postfiltering output
% beam=zeros(M,1);   % beamformer output
ref=zeros(M,2);   % ref output (just two channels, for evaluation of the algorithm
Nref=Nch-1;    % number of reference channels
b=hanning(2*w+1);
b=b/sum(b);
b_xi_local=hanning(2*w_xi_local+1);
b_xi_local=b_xi_local/sum(b_xi_local);
b_xi_global=hanning(2*w_xi_global+1);
b_xi_global=b_xi_global/sum(b_xi_global);
l_mod_lswitch=0;
k_u=round(f_u/Fs*M+1);  % Upper frequency bin for global decision
k_l=round(f_l/Fs*M+1);  % Lower frequency bin for global decision
k_u=min(k_u,M21);
psi0=0.25; 
k_omega=round(250/Fs*M+1);
omega_low=1; omega_high=3;
% omega_low=ones(M21,1); omega_high=repmat(3,M21,1);
%omega_low(1:k_omega)=2; omega_high(1:k_omega)=4;
% omega_low(1:5)=[3 3 8 20 30]; omega_low(6:12)=50; omega_low(13:25)=30;
% omega_low(26:40)=20; omega_low(41:50)=10; omega_low(51:110)=4; omega_low(111:140)=2;
% omega_high=3*omega_low;

r0=1.67; r02=1.81; 
kpsi_u=round(4000/Fs*M+1);  % Upper frequency bin for global decision
kpsi_l=round(1800/Fs*M+1);  % Lower frequency bin for global decision
eta_2termb=1; eta_2term=1; q=0; gamma=1; G=ones(M21,1); xi=0; eta=zeros(M21,1); xi_frame=0; eta_2termr=1; eta_2termz=1;

mu_G = 0.05;
eta_beam = 0.05;
G_beam = zeros(M21,Nref);  
H=ones(M21,Nch);    % blocking matrix
Hnorm = repmat(Nch,M21,1);  % Hnorm = sum(abs(H).^2,2);
Ncount=round(Nframes/10);
fz_idx=cell(Nfin,1);
z=zeros(M,Nch);
z0=zeros(Mno,Nch);

meu=0.1;
% dist=0.05; % inter distance between microphones=5cm
dist=0.1; % inter distance between microphones=5cm
saf_angle=10*pi/180; %degrees
saf1=sin(saf_angle)*(2*pi*dist*Fs)/(330*M);   % sound velocity is 330 m/s
first_time=1; switch_flag=0;
idx8=1:129;
for l=1:Nframes
    
    % Update a frame of the input data
    %--------------------------------------------------------------
    if first_time 
        for k=1:Nfin
            ch_idx=Nch_fin*(k-1)+(1:Nch_fin);
            [z(:,ch_idx),fs,wmode,fz_idx{k}]=readwav(fin{k},'rf',M,0);
        end
    else
        for k=1:Nfin
            ch_idx=Nch_fin*(k-1)+(1:Nch_fin);
            [z0(:,ch_idx),fs,wmode,fz_idx{k}]=readwav(fz_idx{k},'rf',Mno);
        end
        z=[z(Mno+1:M,:); z0];
    end
    Z=fft(win(:,ones(Nch,1)).*z);
    Z=Z(1:M21,:);
    Za2=abs(Z).^2;
    %--------------------------------------------------------------
    
    % esimate stationary background noise of the input channels
    %---------------------------------------------------------------------------
    % smooth over frequency
    Sfz=conv2(Za2,b,'same');  
    if first_time       % first frame
        Sz=Sfz;         
        lambda_davz=Za2;
        Stz=Sfz;
    end
    if l<10       % first few frames
        Sminz=Sz;
        SMactz=Sz;
        SWz=repmat(Sz,[1,1,Nwin]);
        Smintz=Stz;
        Slocaltz=Stz;
        SMacttz=Stz;
        SWtz=repmat(Stz,[1,1,Nwin]);
    end          
    % smooth over time
    Sz=alpha_s*Sz+(1-alpha_s)*Sfz;   
    % Local Minima Search
    Sminz=min(Sminz,Sz);
    SMactz=min(SMactz,Sz);
    I_fz=(Za2<delta_y*Bmin.*Sminz & Sz<delta_s*Bmin.*Sminz);
    conv_Iz=conv2(I_fz,b,'same');
    Sftz=Stz;
    idx=find(conv_Iz);
    if ~isempty(idx)
        if w
            conv_Z=conv2(I_fz.*Za2,b,'same');
            Sftz(idx)=conv_Z(idx)./conv_Iz(idx);
        else
            Sftz(idx)=Za2(idx);
        end
    end
    Stz=alpha_s*Stz+(1-alpha_s)*Sftz;
    Smintz=min(Smintz,Stz);
    SMacttz=min(SMacttz,Stz);
    
    if first_time
        gammaz=ones(M21,Nch);
    else
        gammaz=Za2./max(lambda_dz,1e-10);
    end
    etaz=alpha_eta*eta_2termz+(1-alpha_eta)*max(gammaz-1,0);
    etaz=max(etaz,eta_min);
    vz=gammaz.*etaz./(1+etaz);
    
    qhat=ones(M21,Nch);
    phat=zeros(M21,Nch);
    gamma_mint=Za2./Bmin./max(Smintz,1e-10);
    zetat=Sz./Bmin./max(Smintz,1e-10);
    idx=find(gamma_mint>1 & gamma_mint<delta_yt & zetat<delta_s);
    qhat(idx)=(delta_yt-gamma_mint(idx))/(delta_yt-1);
    phat(idx)=1./(1+qhat(idx)./(1-qhat(idx)).*(1+etaz(idx)).*exp(-vz(idx)));   
    idx=find(gamma_mint>=delta_yt | zetat>=delta_s);
    qhat(idx)=0;
    phat(idx)=1;
    alpha_dt=alpha_d+(1-alpha_d)*phat;
    lambda_davz=alpha_dt.*lambda_davz+(1-alpha_dt).*Za2;
    lambda_dz=Bbeta*lambda_davz;
    
    Gz=ones(M21,Nch);
    idx=find(vz>5);
    Gz(idx)=etaz(idx)./(1+etaz(idx));
    idx=find(vz<=5 & vz>0);
    Gz(idx)=etaz(idx)./(1+etaz(idx)).*exp(0.5*expint(vz(idx)));
    eta_2termz=Gz.^2.*gammaz;
    
    l_mod_lswitch=l_mod_lswitch+1;
    if l_mod_lswitch==Vwin
        switch_flag=1;
        l_mod_lswitch=0;
    else
        switch_flag=0;
    end
    if switch_flag
        SWz=cat(3,SWz(:,:,2:Nwin),SMactz);
        Sminz=min(SWz,[],3);
        SMactz=Sz;
        SWtz=cat(3,SWtz(:,:,2:Nwin),SMacttz);
        Smintz=min(SWtz,[],3);
        SMacttz=Stz;
    end
    %---------------------------------------------------------------------------
    
    % Compute  phase12 and phase12c
    %------------------------------------------------
    if first_time 
        f12=repmat(Z(:,1),1,Nref).*conj(Z(:,2:Nch)); 
    else
        f12=0.6*f12+0.4*repmat(Z(:,1),1,Nref).*conj(Z(:,2:Nch));
    end; 
    phase12=abs(angle(f12))./([ 1 1:(M21-1)]' * (1:Nref));
    phase12c=abs(angle(f12)+angle(H(:,2:Nch)))./([ 1 1:(M21-1)]' * (1:Nref));
    %------------------------------------------------
    
    % Compute  SNR of the input channels
    %-------------------------------------------------------------
    Z_snr=Za2./max(lambda_dz,1e-10);
    %     Z_snr=Za2./max(lambda_davz,1e-10);
    Z_snr=min( repmat(Z_snr(:,1),1,Nref), Z_snr(:,2:Nch) );
    %-------------------------------------------------------------
    
    % update ATF ratios (for the blocking matrix and fixed beamformer)
    %------------------------------------------------------------------------------------------
    for m=2:Nch
        %         leno=sum( Z_snr(17:M21,m-1)>=3 & phase12(17:M21,m-1)<saf1 );
        leno=sum( Z_snr(17:129,m-1)>=3 & phase12(17:129,m-1)<saf1 );
        if leno>40
            %             idx=find(Z_snr(:,m-1)>=4 & phase12(:,m-1)<3*saf1);
            idx=find(Z_snr(:,m-1)>=3 & phase12(:,m-1)<3*saf1);
            amp_calib=Za2(idx,m)./Za2(idx,1);
            if min(amp_calib)>=0.01 & max(amp_calib)<=100;
                H(idx,m)=0.95*H(idx,m)+0.05*Z(idx,m)./Z(idx,1);
            end
        end
    end
    Hnorm = sum(abs(H).^2,2);
    %------------------------------------------------------------------------------------------
    
    % Beamforming
    %---------------------------------------------------------------
    if first_time
        Yr=zeros(M21,Nref);
        YC_buf=Z(:,1);
    end
    Yr_old = Yr;
    YC_buf_old = YC_buf;
    Y_old = YC_buf_old - sum(Yr_old.*G_beam,2);
    
    Yr = Z(:,2:Nch) - H(:,2:Nch) .* repmat(Z(:,1),1,Nref); 
    YC_buf = sum(Z .* conj(H),2) ./ Hnorm; 
    YAS_buf = sum(Yr.*G_beam,2);
    Y = YC_buf - YAS_buf;
    YC_a2=abs(YC_buf).^2;
    Ya2=abs(Y).^2;
    Ya2r=abs(Yr).^2;
    %---------------------------------------------------------------
    
    % Estimate stationary background noise at the beamformer output and the reference channels
    %-----------------------------------------------------------------------------------------------------------------------
    % smooth over frequency
    Sf=conv(b,Ya2);  % 2.1. smooth over frequency
    Sf=Sf(w+1:M21+w); 		 
    if first_time       % first frame
        S=Sf;         
        St=Sf;
    end
    if l<10       % first few frames
        Smin=S;
        SMact=S;
        SW=S*ones(1,Nwin);
        Smint=St;
        Slocalt=St;
        SMactt=St;
        SWt=St*ones(1,Nwin);
    end          
    % smooth over time
    S=alpha_s*S+(1-alpha_s)*Sf;   
    % Local Minima Search
    Smin=min(Smin,S);
    SMact=min(SMact,S);
    I_f=(Ya2<delta_y*Bmin.*Smin & S<delta_s*Bmin.*Smin);
    conv_I=conv(b,I_f);
    conv_I=conv_I(w+1:M21+w); 		 
    Sft=St;
    idx=find(conv_I);
    if ~isempty(idx)
        if w
            conv_Y=conv(b,I_f.*Ya2);
            conv_Y=conv_Y(w+1:M21+w); 		 
            Sft(idx)=conv_Y(idx)./conv_I(idx);
        else
            Sft(idx)=Ya2(idx);
        end
    end
    St=alpha_s*St+(1-alpha_s)*Sft;
    Smint=min(Smint,St);
    SMactt=min(SMactt,St);
    
    % smooth over frequency
    Sfr=conv2(Ya2r,b,'same');  % 2.1. smooth over frequency
    if first_time       % first frame
        Sr=Sfr;         
        lambda_davr=Ya2r;
        Str=Sfr;
    end
    if l<10       % first few frames
        Sminr=Sr;
        SMactr=Sr;
        SWr=repmat(Sr,[1,1,Nwin]);
        Smintr=Str;
        Slocaltr=Str;
        SMacttr=Str;
        SWtr=repmat(Str,[1,1,Nwin]);
    end          
    % smooth over time
    Sr=alpha_s*Sr+(1-alpha_s)*Sfr;   
    % Local Minima Search
    Sminr=min(Sminr,Sr);
    SMactr=min(SMactr,Sr);
    I_fr=(Ya2r<delta_y*Bmin.*Sminr & Sr<delta_s*Bmin.*Sminr);
    conv_Ir=conv2(I_fr,b,'same');
    Sftr=Str;
    idx=find(conv_Ir);
    if ~isempty(idx)
        if w
            conv_Yr=conv2(I_fr.*Ya2r,b,'same');
            Sftr(idx)=conv_Yr(idx)./conv_Ir(idx);
        else
            Sftr(idx)=Ya2r(idx);
        end
    end
    Str=alpha_s*Str+(1-alpha_s)*Sftr;
    Smintr=min(Smintr,Str);
    SMacttr=min(SMacttr,Str);
    
    if first_time
        gamma=ones(M21,1);
    else
        gamma=Ya2./max(Bbeta*lambda_davb,1e-10);
    end
    eta=alpha_eta*eta_2termb+(1-alpha_eta)*max(gamma-1,0);
    eta=max(eta,eta_min);
    v=gamma.*eta./(1+eta);
    
    qhat=ones(M21,1);
    phat=zeros(M21,1);
    gamma_mintY=Ya2./Bmin./max(Smint,1e-10);
    zetat=S./Bmin./max(Smint,1e-10);
    idx=find(gamma_mintY>1 & gamma_mintY<delta_yt & zetat<delta_s);
    qhat(idx)=(delta_yt-gamma_mintY(idx))/(delta_yt-1);
    phat(idx)=1./(1+qhat(idx)./(1-qhat(idx)).*(1+eta(idx)).*exp(-v(idx)));   
    idx=find(gamma_mintY>=delta_yt | zetat>=delta_s);
    qhat(idx)=0;
    phat(idx)=1;
    alpha_dt=alpha_d+(1-alpha_d)*phat;
    if first_time
        lambda_davb=Ya2;
    end
    lambda_davb=alpha_dt.*lambda_davb+(1-alpha_dt).*Ya2;
    
    GH1=ones(M21,1);
    idx=find(v>5);
    GH1(idx)=eta(idx)./(1+eta(idx));
    idx=find(v<=5 & v>0);
    GH1(idx)=eta(idx)./(1+eta(idx)).*exp(0.5*expint(v(idx)));
    eta_2termb=GH1.^2.*gamma;
    
    if first_time
        gammar=ones(M21,Nref);
    else
        gammar=Ya2r./max(lambda_dr,1e-10);
    end
    etar=alpha_eta*eta_2termr+(1-alpha_eta)*max(gammar-1,0);
    etar=max(etar,eta_min);
    vr=gammar.*etar./(1+etar);
    
    qhat=ones(M21,Nref);
    phat=zeros(M21,Nref);
    gamma_mint=Ya2r./Bmin./max(Smintr,1e-10);
    zetat=Sr./Bmin./max(Smintr,1e-10);
    idx=find(gamma_mint>1 & gamma_mint<delta_yt & zetat<delta_s);
    qhat(idx)=(delta_yt-gamma_mint(idx))/(delta_yt-1);
    phat(idx)=1./(1+qhat(idx)./(1-qhat(idx)).*(1+etar(idx)).*exp(-vr(idx)));   
    idx=find(gamma_mint>=delta_yt | zetat>=delta_s);
    qhat(idx)=0;
    phat(idx)=1;
    alpha_dt=alpha_d+(1-alpha_d)*phat;
    if first_time
        lambda_davr=Ya2r;
    end
    lambda_davr=alpha_dt.*lambda_davr+(1-alpha_dt).*Ya2r;
    lambda_dr=Bbeta*lambda_davr;
    
    Gr=ones(M21,Nref);
    idx=find(vr>5);
    Gr(idx)=etar(idx)./(1+etar(idx));
    idx=find(vr<=5 & vr>0);
    Gr(idx)=etar(idx)./(1+etar(idx)).*exp(0.5*expint(vr(idx)));
    eta_2termr=Gr.^2.*gammar;
    
    if switch_flag
        SW=[SW(:,2:Nwin) SMact];
        Smin=min(SW,[],2);
        SMact=S;
        SWt=[SWt(:,2:Nwin) SMactt];
        Smint=min(SWt,[],2);
        SMactt=St;
        
        SWr=cat(3,SWr(:,:,2:Nwin),SMactr);
        Sminr=min(SWr,[],3);
        SMactr=Sr;
        SWtr=cat(3,SWtr(:,:,2:Nwin),SMacttr);
        Smintr=min(SWtr,[],3);
        SMacttr=Str;
    end
    %-----------------------------------------------------------------------------------------------------------------------
    
    % calculate ratio2, ratio2_avr, ratio, ratio1
    %--------------------------------------------------
    %     leniip2=sum(Z_snr>=2)+1;
    %     ratio2=sum(Z_snr>=2 & phase12<saf1)./leniip2;
    leniip2=sum(Z_snr(idx8,:)>=2)+1;
    ratio2=sum(Z_snr(idx8,:)>=2 & phase12(idx8,:)<saf1)./leniip2;
    if first_time, ratio2_avr=ratio2; end
    ratio2_avr=max( ratio2 , 0.5*ratio2_avr+0.5*ratio2 );
    
    %     ratio=sum( Z_snr>=2 & repmat(Ya2<0.25*YC_a2,1,Nref) )./leniip2;
    ratio=sum( Z_snr(idx8,:)>=2 & repmat(Ya2(idx8)<0.25*YC_a2(idx8),1,Nref) )./leniip2;
    if first_time, ratio1=ratio; end
    ratio1=max( ratio, 0.9*ratio1+0.1*ratio );
    %--------------------------------------------------
    
    % update the beam_avr and ref_avr estimators
    %-------------------------------------------------------- 
    if first_time
        beam_avr=repmat(YC_a2,1,Nref);
        ref_avr=beam_avr;
        ref_avr(2:15,:)=ref_avr(2:15,:)*3;
    end;
    for m=1:Nref
        if ratio2_avr(m)<0.5
            idx=find(Z_snr(:,m)>=2 & (phase12(:,m)>saf1*2 | YC_a2<Ya2r(:,m)));
        else
            %             idx=find(Z_snr(:,m)>=4 & phase12(:,m)>saf1 & YC_a2<Ya2r(:,m) & Ya2<=YC_a2);
            idx=find(Z_snr(:,m)>=3 & phase12(:,m)>saf1 & YC_a2<Ya2r(:,m) & Ya2<=YC_a2);
        end
        if length(idx)
            beam_avr(idx,m)=0.9*beam_avr(idx,m)+0.1*YC_a2(idx);
            ref_avr(idx,m)=0.9*ref_avr(idx,m)+0.1*Ya2r(idx,m);
        end
    end
    %-------------------------------------------------------- 
    
    %         % Update G_beam in case only stationry noise is present
    %         if first_time
    %             Pref = sum(Ya2r,2);
    %         else
    %             Pref = lambda_Pref * Pref + (1-lambda_Pref) * sum(Ya2r,2) + eps; 
    %         end
    %         idx=find(r1<r0 & r2<r02);
    %         if length(idx) > M *0.25   % Less than 50% frequency bins contain transients (M/2*0.5)
    %             G_beam(idx,:) = G_beam(idx,:) + mu_G * repmat(Y(idx)./Pref(idx),1,Nref) .* conj(Yr(idx,:));
    %         end
    
    % Update G_beam
    %---------------------------------------------
    if first_time 
        Ya_avr=abs(Y);
    else
        Ya_avr=0.95*Ya_avr+0.05*abs(Y);
    end;
    for m=1:Nref
        if ratio2(m)<0.5
            idx=find( gamma_mintY>1 & ( phase12(:,m)>saf1 | YC_a2<Ya2r(:,m) ) );
        else
            %             idx=find( phase12(:,m)>saf1 & gamma_mintY>4 & YC_a2< Ya2r(:,m) & Ya2<=YC_a2 );
            idx=find( phase12(:,m)>saf1 & gamma_mintY>3 & YC_a2< Ya2r(:,m) & Ya2<=YC_a2 );
        end;
        if length(idx)
            rel=beam_avr(idx,m)./ref_avr(idx,m);
            W1_grad=0.5*rel.*meu.*(conj(Yr_old(idx,m)).*( sign(real(Y_old(idx)))+j*sign(imag(Y_old(idx))))+conj(Yr(idx,m)).*(sign(real(Y(idx)))+j*sign(imag(Y(idx)) ) ) );
            G_beam(idx,m)=G_beam(idx,m)+W1_grad./Ya_avr(idx);
        end;
    end
    %---------------------------------------------
    
    %     % Compute the probability of speech presence, and update the noise estimate at the beamformer output
    %     psi=zeros(M21,1);
    %     qhat=ones(M21,1);
    %     phat=zeros(M21,1);
    %     idx=find(r1>r0 & r2<=r02);
    %     if length(idx)
    %         psi(idx)=1;
    %         qhat(idx)=0;
    %         phat(idx)=1;
    %     end
    %     idx=find(r1>r0 & r2>r02);
    %     if length(idx)
    %         Omega=(S(idx)-Bbeta*lambda_davb(idx))./max(max(Sr(idx,:)-lambda_dr(idx,:),[],2),1e-10);
    %         psi(idx)=max(min((Omega-omega_low(idx))./(omega_high(idx)-omega_low(idx)),1),0);
    %     end
    %     psi_t=mean(psi(kpsi_l:kpsi_u));
    %     
    %     if psi_t<=psi0  % in this case the frame does not contain speech
    %         qhat=ones(M21,1);
    %         phat=zeros(M21,1);
    %     elseif length(idx)   % in this case the frame might contain speech
    %         qhat(idx)=min(max(max((delta_yt-gamma_mint(idx))/(delta_yt-1),1-psi(idx)),0),0.999);
    %         phat(idx)=1./(1+qhat(idx)./(1-qhat(idx)).*(1+eta(idx)).*exp(-v(idx)));   
    %     end
    %     alpha_dt=alpha_d+(1-alpha_d)*phat;
    %     lambda_dav=alpha_dt.*lambda_dav+(1-alpha_dt).*Ya2;
    %     % 2.4. Noise Spectrum Estimate
    %     %lambda_d=Bbeta*lambda_dav;
    %     if Nseg_k>8
    %         lambda_d=noise_factor*lambda_dav;
    %     else
    %         lambda_d=Bbeta*lambda_davb;
    %     end
    
    % Compute  the LNS (local non-stationarity) of the beamformer output and the reference channels
    %----------------------------------------------------------------------------------------------------------------------------
    r1=S/Bbeta./max(lambda_davb,1e-10);
    r2=max(max(Sr./max(lambda_dr,1e-10),[],2),1e-10);
    %----------------------------------------------------------------------------------------------------------------------------
    
    % Update the nonstationary noise estimate at the beamformer output
    %-----------------------------------------------------------------------------------------------
    if first_time
        gamma=ones(M21,1);
    else
        gamma=Ya2./max(lambda_d,1e-10);
    end
    eta=alpha_eta*eta_2term+(1-alpha_eta)*max(gamma-1,0);
    eta=max(eta,eta_min);
    v=gamma.*eta./(1+eta);
    
    %     % This section is for aggressive noise reduction
    %     %----------------------------------------------------------------
    %     psi=zeros(M21,1);
    %     qhat=ones(M21,1);
    %     phat=zeros(M21,1);
    %     H1_flag= ~sum(Z_snr>=3 &  phase12>3*saf1,2);
    %     idx=find(r1>r0 & r2<=r02 & H1_flag);
    %     %     idx=find(r1>r0 & r2<=r02);
    %     if length(idx)
    %         psi(idx)=1;
    %         qhat(idx)=0;
    %         phat(idx)=1;
    %     end
    %     idx=find(r1>r0 & r2>r02 & H1_flag);
    %     %     idx=find(r1>r0 & r2>r02);
    %     if length(idx)
    %         Omega=(S(idx)-Bbeta*lambda_davb(idx))./max(max(Sr(idx,:)-lambda_dr(idx,:),[],2),1e-10);
    %         psi(idx)=max(min((Omega-omega_low)./(omega_high-omega_low),1),0);
    %         qhat(idx)=min(max(max((delta_yt-gamma_mintY(idx))/(delta_yt-1),1-psi(idx)),0),0.999);
    %         phat(idx)=1./(1+qhat(idx)./(1-qhat(idx)).*(1+eta(idx)).*exp(-v(idx)));   
    %     end
    %     %----------------------------------------------------------------
    
    % This section is for non-aggressive noise reduction
    %----------------------------------------------------------------
    qhat=zeros(M21,1);
    phat=ones(M21,1);
    idx=find(r1<=r0);
    if length(idx)
        qhat(idx)=1;
        phat(idx)=0;
    end
    idx=find(r1>r0 & r2>r02 & phase12>saf1);
    if length(idx)
        Omega=(S(idx)-Bbeta*lambda_davb(idx))./max(max(Sr(idx,:)-lambda_dr(idx,:),[],2),1e-10);
        psi=max(min((Omega-omega_low)./(omega_high-omega_low),1),0);
        qhat(idx)=min(max(max((delta_yt-gamma_mintY(idx))/(delta_yt-1),1-psi),0),0.999);
        phat(idx)=1./(1+qhat(idx)./(1-qhat(idx)).*(1+eta(idx)).*exp(-v(idx)));   
    end
    %----------------------------------------------------------------
    
    alpha_dt=alpha_d+(1-alpha_d)*phat;
    if first_time
        lambda_dav=Ya2;
    end
    lambda_dav=alpha_dt.*lambda_dav+(1-alpha_dt).*Ya2;
    % Noise Spectrum Estimate
    %        lambda_d=Bbeta*lambda_davb;
    lambda_d=noise_factor*lambda_dav;
    
    %     if first_time, lambda_d=lambda_davb; end; 
    %     index=find(Ya2<5*lambda_davb & phase12>saf1);
    %     lambda_d(index)=0.95*lambda_d(index)+0.05*Ya2(index);
    
    %-----------------------------------------------------------------------------------------------
    
    % Estimate the a priori probability for signal absence
    %--------------------------------------------------------------------------
    xi=alpha_xi*xi+(1-alpha_xi)*eta;
    xi_local=conv(xi,b_xi_local);
    xi_local=xi_local(w_xi_local+1:M21+w_xi_local);
    xi_global=conv(xi,b_xi_global);
    xi_global=xi_global(w_xi_global+1:M21+w_xi_global);
    dxi_frame=xi_frame;
    xi_frame=mean(xi(k_l:k_u));
    dxi_frame=xi_frame-dxi_frame;
    if xi_local>0, xi_local_dB=10*log10(xi_local); else, xi_local_dB=-100; end
    if xi_global>0, xi_global_dB=10*log10(xi_global); else, xi_global_dB=-100; end
    if xi_frame>0, xi_frame_dB=10*log10(xi_frame); else, xi_frame_dB=-100; end
    
    P_local=ones(M21,1);
    idx=find(xi_local_dB<=xi_ll_dB);
    P_local(idx)=P_min;
    idx=find(xi_local_dB>xi_ll_dB & xi_local_dB<xi_lu_dB);
    P_local(idx)=P_min+(xi_local_dB(idx)-xi_ll_dB)/(xi_lu_dB-xi_ll_dB)*(1-P_min);
    
    P_global=ones(M21,1);
    idx=find(xi_global_dB<=xi_gl_dB);
    P_global(idx)=P_min;
    idx=find(xi_global_dB>xi_gl_dB & xi_global_dB<xi_gu_dB);
    P_global(idx)=P_min+(xi_global_dB(idx)-xi_gl_dB)/(xi_gu_dB-xi_gl_dB)*(1-P_min);
    
    if xi_frame_dB<=xi_fl_dB
        P_frame=P_min;
    elseif dxi_frame>=0
        xi_m_dB=min(max(xi_frame_dB,xi_ml_dB),xi_mu_dB);
        P_frame=1;
    elseif xi_frame_dB>=xi_m_dB+xi_fu_dB
        P_frame=1;
    elseif xi_frame_dB<=xi_m_dB+xi_fl_dB
        P_frame=P_min;
    else
        P_frame=P_min+(xi_frame_dB-xi_m_dB-xi_fl_dB)/(xi_fu_dB-xi_fl_dB)*(1-P_min);
    end   
    
    q=1-P_global.*P_local*P_frame;
    q=min(q,q_max);
    %--------------------------------------------------------------------------
    
    % Compute the spectral gain
    %--------------------------------------------------------------------------
    gamma=Ya2./max(lambda_d,1e-10);
    eta=alpha_eta*eta_2term+(1-alpha_eta)*max(gamma-1,0);
    eta=max(eta,eta_min);
    v=gamma.*eta./(1+eta);
    % Conditional Probability for Signal-Presence Estimate
    %     PH1=1./(1+q./(1-q).*(1+eta).*exp(-v));
    PH1=zeros(M21,1);
    idx=find(q<0.9);
    PH1(idx)=1./(1+q(idx)./(1-q(idx)).*(1+eta(idx)).*exp(-v(idx)));
    
    % 7. Spectral Gain
    GH1=ones(M21,1);
    idx=find(v>5);
    GH1(idx)=eta(idx)./(1+eta(idx));
    idx=find(v<=5 & v>0);
    GH1(idx)=eta(idx)./(1+eta(idx)).*exp(0.5*expint(v(idx)));
    eta_2term=GH1.^2.*gamma;
    
    %     G=GH1.^PH1.*G_f.^(1-PH1);
    %     G=max(G,G_f);
    
    Smint_global=[Smint [Smint(2:M21);Smint(M21)] [Smint(3:M21);Smint(M21-1:M21)] [Smint(4:M21);Smint(M21-2:M21)] [Smint(1);Smint(1:M21-1)] [Smint(1:2);Smint(1:M21-2)] [Smint(1:3);Smint(1:M21-3)]];
    Smint_global=min(Smint_global,[],2);
    lambda_d_global=1.5*Bmin*Smint_global;
    if first_time
        Sy=Ya2;
    else
        Sy=0.8*Sy+0.2*Ya2; 
    end
    GH0=G_f*(lambda_d_global./Sy).^0.5;
    %     GH0=G_f*(Bbeta*lambda_davb./Sy).^0.5;
    
    G=GH1.^PH1.*GH0.^(1-PH1);
    %--------------------------------------------------------------------------
    
    % Modify the gain
    %--------------------------------------------------
    %     leno2=sum( gamma_mintY>2  & Ya2>0.64*YC_a2 )+1;
    %     leno1_ins=sum( repmat(gamma_mintY>2  & Ya2>0.64*YC_a2,1,Nref) & phase12<saf1 );
    leno2=sum( gamma_mintY(idx8)>2  & Ya2(idx8)>0.64*YC_a2(idx8) )+1;
    leno1_ins=sum( repmat(gamma_mintY(idx8)>2  & Ya2(idx8)>0.64*YC_a2(idx8),1,Nref) & phase12(idx8,:)<saf1 );
    if first_time 
        leno1=leno1_ins;
    end;
    leno1=0.9*leno1+0.1*leno1_ins;
    ratio3=leno1/leno2;
    if first_time
        ratio3_avr=ratio3;
    end;
    ratio3_avr=0.95*ratio3_avr+0.05*ratio3;
    for m=1:Nref
        if ratio1(m)<0.15 & ratio2_avr(m)>0.5
            idx=find(Z_snr(:,m)>=2 & G>0.5 & Ya2>0.64*YC_a2);
            %             idx=find(Z_snr(:,m)>=3 & G>0.5 & Ya2>0.7*YC_a2);
            %             G(idx)=max(G(idx),1);
            G(idx)=GH1(idx);
            %                         G(idx)=max(GH1(idx),1);
        end;
        if ratio1(m)>0.15
            idx=find(Ya2<0.25*YC_a2);
            %             G(idx)=G_f;
            G(idx)=GH0(idx);
        end;
        if (ratio2_avr(m)<0.5 & ratio3_avr(m)<0.4 & leno1<90)  
            %             idx=find(eta_2term>G_f^2);
            %             G(idx)=G_f./gamma(idx).^0.5;
            G=GH0;
            %             idx=find(phase12(:,m)<saf1 & eta_2term>G_f^2);
            %             G(idx)=max(G(idx),G_f);
        else
            idx=17:M21;
            %             pidx=find( (phase12(idx,m)>2*saf1  &   eta_2term(idx)>4) | (phase12(idx,m)>saf1  & eta_2term(idx)>8) );
            pidx=find( (phase12c(idx,m)>2*saf1  &   eta_2term(idx)>3) | (phase12c(idx,m)>saf1  & eta_2term(idx)>5) );
            %             pidx=find( (phase12(idx,m)>2*saf1  &   eta_2term(idx)>3) | (phase12(idx,m)>saf1  & eta_2term(idx)>5) );  % this produces practically similar results as using phase12c
            G(idx(pidx))=G(idx(pidx))/4;
            idx=find( phase12(:,m)<1.5*saf1 & eta_2term>3 );
            %             G(idx)=max(G(idx),1);
            G(idx)=GH1(idx);
        end;
    end
    %--------------------------------------------------
    
    X=[zeros(2,1); G(3:M21-1).*Y(3:M21-1); 0]; 
    X(M21+1:M)=conj(X(M21-1:-1:2)); %extend the anti-symmetric range of the spectum 
    x=Cwin^2*win.*real(ifft(X));
    out=out+x;
    Y(M21+1:M)=conj(Y(M21-1:-1:2)); %extend the anti-symmetric range of the spectum 
    y=Cwin^2*win.*real(ifft(Y));
    %     beam=beam+y;
    if first_time
        if out_st
            fout_idx=writewav(out(1:Mno,[1 1]) ,fs,fout,'rf',0);    % stereo output
        else
            fout_idx=writewav(out(1:Mno),fs,fout,'rf',0);
        end
        %         fbeam_idx=writewav(beam(1:Mno),fs,fbeam,'rf',0);
    else
        if out_st
            fout_idx=writewav(out(1:Mno,[1 1]),fs,fout_idx,'rf');         % stereo output
        else
            fout_idx=writewav(out(1:Mno),fs,fout_idx,'rf');         
        end
        %         fbeam_idx=writewav(beam(1:Mno),fs,fbeam_idx,'rf');
    end
    if exist('fref')
        Yr_tmp=Yr(:,1:2); % take two ref channels, for evaluation
        Yr_tmp(M21+1:M,:)=conj(Yr_tmp(M21-1:-1:2,:));
        yr=Cwin^2*win(:,[1 1]).*real(ifft(Yr_tmp));
        ref=ref+yr;
        if first_time
            fref_idx=writewav(ref(1:Mno,:),fs,fref,'rf',0);
        else
            fref_idx=writewav(ref(1:Mno,:),fs,fref_idx,'rf');
        end
        ref=[ref(Mno+1:M,:); zeros(Mno,2)];
    end
    out=[out(Mno+1:M); zeros(Mno,1)];
    %     beam=[beam(Mno+1:M); zeros(Mno,1)];
    if ~mod(l,Ncount), disp([sprintf('%3.0f',100*l/Nframes) '%']); end
    first_time=0;
end
if out_st
    fout_idx=writewav(out(1:M-Mno,[1 1]),fs,fout_idx,'rf',0);        % stereo output
else
    fout_idx=writewav(out(1:M-Mno),fs,fout_idx,'rf',0);
end
% fbeam_idx=writewav(beam(1:M-Mno),fs,fbeam_idx,'rf',0);
fclose(fout_idx(1));
% fclose(fbeam_idx(1));
for k=1:Nfin
    fclose(fz_idx{k}(1));
end
if exist('fref')
    fref_idx=writewav(ref(1:M-Mno,:),fs,fref_idx,'rf',0);
    fclose(fref_idx(1));
end

function y = lnshift(x,t)
% lnshift -- t circular left shift of 1-d signal
%  Usage
%    y = lnshift(x,t)
%  Inputs
%    x   1-d signal
%  Outputs
%    y   1-d signal 
%        y(i) = x(i+t) for i+t < n
%	 		y(i) = x(i+t-n) else
% ***************************************************************@
szX=size(x);
if szX(1)>1
    n=szX(1);
    y=[x((1+t):n); x(1:t)];
else
    n=szX(2);
    y=[x((1+t):n) x(1:t)];
end