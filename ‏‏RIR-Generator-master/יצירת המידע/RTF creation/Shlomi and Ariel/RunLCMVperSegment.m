function [yO,yOO,w,SPP] = RunLCMVperSegment(basePath,fs,mref,z,seat1_16kHz,seat3_16kHz,seat4_16kHz,vadSpk,vadInt4,Tstat,vadSeat1,vadSeat3,vadSeat4,resultDir,SPP)

[M,Lz]  = size(z);
% M:    number of microphones
% Lz:   length of noisy signal

%% LCMV parameters:
g               = [1 ; 0];             % (2 x 1)    desired response vector
alpha_R         = 0.95;                % (1 x 1)    Rvv update parameter
rho_th          = 150;                 % (1 x 1)    threshold of SPP decision
sum_th          = 257;
NBuffer         = 16;                  % (1 x 1)    size of segments for update RTF
segCounter      = 0;                   % (1 x 1)    counter of segments

%% mix-max parameters:
addpath([basePath,'voicebox']);
load([basePath,'MixMaxInputs.mat']);
mu_i_k           = MixMax.mu_i_k;      % (K x M)  GMM means (mus)
K                = size(mu_i_k,1);     % (1 x 1)  number of samples in frequency domain (K = Nfft/2 + 1)
NumMixtures      = size(mu_i_k,2);     % (1 x 1)  number of mixtures
Nfft             = 2*(K-1);            % (1 x 1)
sigma_i_k        = MixMax.sigma_i_k;   % (K x M)  GMM standard deviations (sigmas)
synt_win         = MixMax.synt_win.';  % (L x 1)  synthesis window
b1               = MixMax.b1;
b2               = MixMax.b2;
W1               = MixMax.W1;
W2               = MixMax.W2;
beta             = 0.06;               % (1 x 1)  noise reduction level (eq. 16)
alpha            = 0.004*ones(K,1);    % (K x 1)  smoothing parameter of the noise update (eq. 25)
ovrlp            = 0.75;               % (1 x 1)  overlap between following frames
sub_num          = 1/(1 - ovrlp) - 1;  % (1 x 1)
n1               = 8;                  % (1 x 1)  number of neighbors
analysis_win     = [hanning(Nfft-1); 0].';
noise_adaptation = 1;                  %  flag    update noise

%% Mix-Max - Initializations:
v           = z(:,Tstat); % (Lv x 1) noise segment
Lv       	= size(v,2);
Nseg        = fix(Lz/(Nfft*(1 - ovrlp))) - sub_num; % number of segments of the noisy signal
noiseNseg	= fix(Lv/(Nfft*(1 - ovrlp))) - sub_num; % number of segments of the noise signal

% estimate the gain of the speech signal:
xGain = std(z(mref,:))^2 - std(v(mref,:))^2;
if xGain < 0
    xGain = std(z(mref,:))*0.5;
else
    xGain = sqrt(xGain);
end
logxGain = log(xGain);

% find the probabilities of each Gaussian for all the signal and scale
% each vector to the range [0-1]: C = (C - mean(C))/std(C))
C = melcepst(z(mref,:)/std(z(mref,:)),fs,'e0dD',12,floor(3*log(fs)),Nfft,Nfft*(1-ovrlp),0,0.5);
C = (C - repmat(mean(C),size(C,1),1))./(repmat(std(C),size(C,1),1));

% add context frames to preserve smoothnes
sizeC = size(C);
n = n1*2 + 1;
temp_mat = zeros(sizeC(1) + n1*2,sizeC(2));
temp_mat(n1+1:size(temp_mat,1)-n1,:) = C;
temp_C = zeros(sizeC(1),sizeC(2)*n);
for i = (n1 + 1):(size(temp_mat,1) - n1)
    for j = 1:n
        t = (1:sizeC(2)) + (sizeC(2)*(j-1));
        temp_C(i - n1,t) = temp_mat(i - n1 - 1 + j,:);
    end
end
C = temp_C';

% feed-forward in the NN (updated in 19/06/2016):
A = W1*C + repmat(b1,1,size(C,2));
H1 = max(A,0)*0.5;
B = W2*H1 + repmat(b2,1,size(C,2));
P_mat = (exp(B))./repmat(sum(exp(B),1),39,1);

% Noise spectrum. find the STFT of the first samples (speech free)
logY = zeros(K,noiseNseg); % (K x noiseNseg)   log-spectrum of the noise signal
for segIdx = 1:noiseNseg
    TimeVec = (segIdx-1)*Nfft*(1-ovrlp)+1 : (segIdx-1)*Nfft*(1-ovrlp)+Nfft; % (1 x Nfft). frame overlapping of 75%
    V = fft(v(mref,TimeVec).*analysis_win,Nfft,2); % (Nfft x 1). spectrum of noise segment
    logY(:,segIdx) = log(abs(V(1:K))); % (Nfft x 1). log-spectrum of the noise signal
end

% find noise's parameters
mu_y = mean(logY,2);  % (K x 1). noise means (for each frequency sample)
sigma_y = std(logY,0,2); % (K x 1). noise standard deviations (for each frequency sample)

logfConst = -log(2*pi)/2 - log(sigma_i_k); % (K x NumMixtures). eq. (4): constant part of log[f_i_k(X_k)]

%% LCMV - Initializations:

% calc initial staionary noise PSD (Rvv)
Rvv = calcR(v,Nfft,Nfft/2);
eyeMat = repmat(eye(M),[1,1,Nfft]);

% calc seat1 RTF
seat1Mat = seat1_16kHz(:,vadSeat1 == 1);
Rs1s1 = calcR(seat1Mat,Nfft,Nfft/2); % calc seat1 PSD
vv = gevd_freq(Rs1s1,eyeMat);
RTFseat1 = squeeze(vv(:,1,:)./repmat(vv(mref,1,:),[M,1,1]));

% calc seat3 RTF
% seat3Mat = seat3_16kHz(:,vadSeat3 == 1);
% Rs3s3 = calcR(seat3Mat,Nfft,Nfft/2); % calc seat3 PSD
% vv = gevd_freq(Rs3s3,eyeMat);
% RTFseat3 = squeeze(vv(:,1,:)./repmat(vv(mref,1,:),[M,1,1]));

% calc seat4 RTF
seat4Mat = seat4_16kHz(:,vadSeat4 == 1);
Rs4s4 = calcR(seat4Mat,Nfft,Nfft/2); % calc seat4 PSD
vv = gevd_freq(Rs4s4,eyeMat);
RTFseat4 = squeeze(vv(:,1,:)./repmat(vv(mref,1,:),[M,1,1]));

% calc initial RTF of desired speaker
% speakerMat = z(:,vadSpk == 1);
% Rss = calcR(speakerMat,Nfft,Nfft/2); % calc desired speaker PSD
% vv = gevd_freq(Rss,Rvv);
% RTFdesd = squeeze(vv(:,1,:)./repmat(vv(mref,1,:),[M,1,1]));
RTFdesd = RTFseat1;

% calc initial RTF of interference3
% interferenceMat = z(:,vadInt3 == 1);
% Ri3i3 = calcR(interferenceMat,Nfft,Nfft/2); % calc interference3 PSD
% vv = gevd_freq(Ri3i3,Rvv);
% RTFint3 = squeeze(vv(:,1,:)./repmat(vv(mref,1,:),[M,1,1]));

% calc initial RTF of interference4
% interferenceMat = z(:,vadInt4 == 1);
% Ri4i4 = calcR(interferenceMat,Nfft,Nfft/2); % calc interference4 PSD
% vv = gevd_freq(Ri4i4,Rvv);
% RTFint4 = squeeze(vv(:,1,:)./repmat(vv(mref,1,:),[M,1,1]));
RTFint4 = RTFseat4;

% calc initial RTF of double talk
% doubleTalkMat = z(:,9*fs+1:12*fs);
% Rdtdt = calcR(doubleTalkMat,Nfft,Nfft/2); % calc double talk PSD
% vv = gevd_freq(Rdtdt,Rvv);
% RTFdt = squeeze(vv(:,1,:)./repmat(vv(mref,1,:),[M,1,1]));

% debug: check cos_teta RTFs
% --------------------------
% cos_teta1_desd = zeros(Nfft,1);
% cos_teta3_desd = zeros(Nfft,1);
% cos_teta4_desd = zeros(Nfft,1);
% for f = 1:Nfft
%     cos_teta1_desd(f) = abs(dot(RTFdesd(:,f),RTFseat1(:,f))/(norm(RTFdesd(:,f))*norm(RTFseat1(:,f))));
%     cos_teta3_desd(f) = abs(dot(RTFdesd(:,f),RTFseat3(:,f))/(norm(RTFdesd(:,f))*norm(RTFseat3(:,f))));
%     cos_teta4_desd(f) = abs(dot(RTFdesd(:,f),RTFseat4(:,f))/(norm(RTFdesd(:,f))*norm(RTFseat4(:,f))));
% end
% 
% h = figure; hold on;
% plot(cos_teta1_desd,'r');
% plot(cos_teta3_desd,'g');
% plot(cos_teta4_desd,'b');
% title('cos\_teta: desired speaker');
% xlabel('frequency bin');
% legend('seat1','seat3','seat4');
% axis([0 length(cos_teta1_desd) 0 1]);
% hold off;
% savefig([resultDir,'\cos_teta_desd.fig']);
% close(h);
% 
% [~,idx] = max([cos_teta1_desd,cos_teta3_desd,cos_teta4_desd],[],2);
% aaa1 = idx;
% aaa1(aaa1 ~= 1) = nan;
% aaa2 = idx;
% aaa2(aaa2 ~= 2) = nan;
% aaa3 = idx;
% aaa3(aaa3 ~= 3) = nan;
% 
% h = figure; hold on;
% stem(aaa1,'r');
% stem(aaa2,'g');
% stem(aaa3,'b');
% title('max(cos\_teta): desired speaker');
% xlabel('frequency bin');
% legend('seat1','seat3','seat4');
% axis([0 length(cos_teta1_desd) 0 4]);
% savefig([resultDir,'\cos_teta_desd_max.fig']);
% close(h);
% 
% % sum(idx == 1)
% % sum(idx == 2)
% % sum(idx == 3)
% 
% cos_teta1_Int4 = zeros(Nfft,1);
% cos_teta3_Int4 = zeros(Nfft,1);
% cos_teta4_Int4 = zeros(Nfft,1);
% for f = 1:Nfft
%     cos_teta1_Int4(f) = abs(dot(RTFint4(:,f),RTFseat1(:,f))/(norm(RTFint4(:,f))*norm(RTFseat1(:,f))));
%     cos_teta3_Int4(f) = abs(dot(RTFint4(:,f),RTFseat3(:,f))/(norm(RTFint4(:,f))*norm(RTFseat3(:,f))));
%     cos_teta4_Int4(f) = abs(dot(RTFint4(:,f),RTFseat4(:,f))/(norm(RTFint4(:,f))*norm(RTFseat4(:,f))));
% end
% 
% h = figure; hold on;
% plot(cos_teta1_Int4,'r');
% plot(cos_teta3_Int4,'g');
% plot(cos_teta4_Int4,'b');
% title('cos\_teta: interference 4');
% xlabel('frequency bin');
% legend('seat1','seat3','seat4');
% axis([0 length(cos_teta1_Int4) 0 1]);
% hold off;
% savefig([resultDir,'\cos_teta_int4.fig']);
% close(h);
% 
% [~,idx] = max([cos_teta1_Int4,cos_teta3_Int4,cos_teta4_Int4],[],2);
% aaa1 = idx;
% aaa1(aaa1 ~= 1) = nan;
% aaa2 = idx;
% aaa2(aaa2 ~= 2) = nan;
% aaa3 = idx;
% aaa3(aaa3 ~= 3) = nan;
% 
% h = figure; hold on;
% stem(aaa1,'r');
% stem(aaa2,'g');
% stem(aaa3,'b');
% title('max(cos\_teta): interference 4');
% xlabel('frequency bin');
% legend('seat1','seat3','seat4');
% axis([0 length(cos_teta1_Int4) 0 4]);
% savefig([resultDir,'\cos_teta_Int4_max.fig']);
% close(h);
% 
% % sum(idx == 1)
% % sum(idx == 2)
% % sum(idx == 3)
% 
% cos_teta1_dt = zeros(Nfft,1);
% cos_teta3_dt = zeros(Nfft,1);
% cos_teta4_dt = zeros(Nfft,1);
% for f = 1:Nfft
%     cos_teta1_dt(f) = abs(dot(RTFdt(:,f),RTFseat1(:,f))/(norm(RTFdt(:,f))*norm(RTFseat1(:,f))));
%     cos_teta3_dt(f) = abs(dot(RTFdt(:,f),RTFseat3(:,f))/(norm(RTFdt(:,f))*norm(RTFseat3(:,f))));
%     cos_teta4_dt(f) = abs(dot(RTFdt(:,f),RTFseat4(:,f))/(norm(RTFdt(:,f))*norm(RTFseat4(:,f))));
% end
% 
% h = figure; hold on;
% plot(cos_teta1_dt,'r');
% plot(cos_teta3_dt,'g');
% plot(cos_teta4_dt,'b');
% title('cos\_teta: double talk');
% xlabel('frequency bin');
% legend('seat1','seat3','seat4');
% axis([0 length(cos_teta1_dt) 0 1]);
% hold off;
% savefig([resultDir,'\cos_teta_dt.fig']);
% close(h);
% 
% [~,idx] = max([cos_teta1_dt,cos_teta3_dt,cos_teta4_dt],[],2);
% aaa1 = idx;
% aaa1(aaa1 ~= 1) = nan;
% aaa2 = idx;
% aaa2(aaa2 ~= 2) = nan;
% aaa3 = idx;
% aaa3(aaa3 ~= 3) = nan;
% 
% h = figure; hold on;
% stem(aaa1,'r');
% stem(aaa2,'g');
% stem(aaa3,'b');
% title('max(cos\_teta): Double Talk');
% xlabel('frequency bin');
% legend('seat1','seat3','seat4');
% axis([0 length(cos_teta1_dt) 0 4]);
% savefig([resultDir,'\cos_teta_dt_max.fig']);
% close(h);
% 
% % sum(idx == 1)
% % sum(idx == 2)
% % sum(idx == 3)

%% main loop:
w           = zeros(M,Nfft,Nseg);    % (M x Nfft x Nseg)    beamformer coeffs
yO          = zeros(1,Lz);           % (1 x Lz)             LCMV output
yOO         = zeros(1,Lz);           % (1 x Lz)             Mix-Max output
zbufForRTF  = zeros(M,NBuffer*Nfft); % (M x NBuffer*Nfft)  	buffer for RTF
% if ~exist('SPP','var')
%     SPP     = zeros(K,Nseg);         % (K x Nseg)        speech presence probability (rho)
% end
sumSPP      = zeros(1,Nseg);
for segIdx = 1:Nseg
    
    TimeVec = (segIdx-1)*Nfft*(1-ovrlp)+1 : (segIdx-1)*Nfft*(1-ovrlp) + Nfft; % (1 x Nfft). frame overlapping of 75%
    Zbuf = fft(z(:,TimeVec).*repmat(analysis_win,M,1),Nfft,2); % input segment of z (freq)
    zbufForRTF = [zbufForRTF(:,1:(NBuffer-1)*Nfft) , z(:,TimeVec)];
      
    if (segIdx > 1)
        
        sumSPP(segIdx) = sum(SPP(:,segIdx-1)); % hard decision
        
        if (sumSPP(segIdx) <= rho_th)
            
            % update noise
            for f = 1:Nfft
                Rvv(:,:,f) = alpha_R*(Zbuf*Zbuf') + (1-alpha_R)*Rvv(:,:,f);
            end
            segCounter = 0;
            
        else
            
            % update speaker
            if (segCounter < NBuffer)
                
                segCounter = segCounter + 1;
                
            else
                
                % calc current RTF
                Rcc = calcR(zbufForRTF,Nfft,Nfft/2); % calc current PSD
                vv = gevd_freq(Rcc,Rvv);
                RTFcurrent = squeeze(vv(:,1,:)./repmat(vv(mref,1,:),[M,1,1]));
                
                % calc cosine distance:
                costeta1 = zeros(Nfft,1);
                % costeta3 = zeros(Nfft,1);
                costeta4 = zeros(Nfft,1);
                for f = 1:Nfft
                    costeta1(f) = abs(dot(RTFcurrent(:,f),RTFseat1(:,f))/(norm(RTFcurrent(:,f))*norm(RTFseat1(:,f))));
                    % costeta3(f) = abs(dot(RTFcurrent(:,f),RTFseat3(:,f))/(norm(RTFcurrent(:,f))*norm(RTFseat3(:,f))));
                    costeta4(f) = abs(dot(RTFcurrent(:,f),RTFseat4(:,f))/(norm(RTFcurrent(:,f))*norm(RTFseat4(:,f))));
                end
                [~,idx] = max([costeta1,costeta4],[],2);
                sum1 = sum(idx == 1);
                % sum3 = sum(idx == 2);
                sum4 = sum(idx == 2);
                if (sum1 > sum4) && (sum1 > sum_th)
                    RTFdesd = RTFcurrent;
                % elseif (sum3 > sum1) && (sum3 > sum4) && (sum3 > sum_th)
                %     RTFint3 = RTFcurrent;
                elseif (sum4 > sum1) && (sum4 > sum_th)
                    RTFint4 = RTFcurrent;
                else
                    segCounter = 0;
                end
                
            end % if (segCounter < NBuffer)
            
        end % if (sumSPP(segIdx) < rho_th)
        
    end % if (segIdx > 1)
    
    for f = 1:Nfft
        
        % calc w
        C1 = [RTFdesd(:,f),RTFint4(:,f)];
        w(:,f,segIdx) = (Rvv(:,:,f)^-1)*C1*((C1'*(Rvv(:,:,f)^-1)*C1)^-1)*g;
        
    end
    
    % apply beamformer:
    % -----------------
    YO = dot(w(:,:,segIdx),Zbuf,1);
    yObuf = real(ifft(YO,Nfft).*synt_win); % back to time domain
    yO(TimeVec) = yO(TimeVec) + yObuf(1:Nfft); % insert segment to output (overlap & add)
    
    % apply mix-max:
    % --------------
    AmpYO = log(abs(YO(1:K).')); % (K x 1). log(spectrum Amplitude)
    PhaseYO = angle(YO.'); % (Nfft x 1). phase
    
    logf_exp = (AmpYO*ones(1,NumMixtures) - mu_i_k - logxGain)./sigma_i_k; % (K x NumMixtures). eq. (4): exponent part of log[f_i_k(X_k)]
    logf = logfConst - logf_exp.*logf_exp/2; % (K x NumMixtures). eq. (4): log[f_i_k(X_k)]
    logF = log(0.5 + 0.5*erf(logf_exp/sqrt(2)) + eps); % (K x NumMixtures). Cumulative Distribution Function of f_i_k(X_k)
    
    logg_exp = (AmpYO - mu_y) ./ sigma_y; % (K x 1). eq. (6): exponent part of log[g_k(Y_k)]
    logg = -log(2*pi)/2 - log(sigma_y) - logg_exp.*logg_exp/2; % (K x 1). eq. (6): log[g_k(Y_k)]
    logG = log(0.5 + 0.5*erf(logg_exp/sqrt(2)) + eps); % (K x 1). cumulative distribution function of log[g_k(Y_k)]
    
    R_i_k = exp(logf - logF); % (K x NumMixtures). eq. (12)
    R_Y_k = exp(logg - logG); % (K x 1). eq. (12)
    rho_i_k = R_i_k ./ (R_i_k + R_Y_k*ones(1,NumMixtures) + eps); % (K x NumMixtures). eq. (12)
    
    % Calculate the Speech Presense Probability (SPP):
    SPP(:,segIdx) = rho_i_k*P_mat(:,segIdx); % (K x 2258). eq. (19)
    
    % The processing:
    AmpYOO = (exp(AmpYO).^(SPP(:,segIdx))).*((beta*exp(AmpYO)).^(1 - SPP(:,segIdx))); % (K x 1). eq. (17)
    
    % Update Noise Estimation:
    if (noise_adaptation == 1) && (segIdx > noiseNseg)
        mu_y = SPP(:,segIdx).*mu_y + (1-SPP(:,segIdx)).*(alpha.*AmpYO + (1 - alpha).*mu_y); % (K x 1). eq. (25)
        sigma_y = SPP(:,segIdx).*sigma_y + (1-SPP(:,segIdx)).*(alpha.*abs(AmpYO - mu_y) + (1 - alpha).*sigma_y); % (K x 1). eq. (25)
    end
    
    % Synthesis:
    YOO = AmpYOO([1:K,K-1:-1:2]).*exp(1i*PhaseYO);
    yOObuf = (real(ifft(YOO,Nfft)).*synt_win.').';
    yOO(TimeVec) = yOO(TimeVec) + yOObuf(1:Nfft); % insert segment to output (overlap & add)
    
end % segment

end % function RunLCMVperSegment