function Rxx = calcR_mul(x,Nfft,JUMP)

% Estimate PSD (Power spectral density) using Welch method

% inputs:
% -------
% x:    noise segment (M microphones)
% M:    number of Microphones
% N:    number of samples
% Lx:   length of x
% Nfft: fft length
% Q:    number of frames
%
% outputs:
% --------
% Rxx:   spectrum of the noise (M x M x Nfft)

%% initializations
[M,Lx,N] = size(x);
SegNo = floor((Lx - (Nfft - JUMP))/JUMP); % number of segments
win = hamming(Nfft);                    % analysis window
X = zeros(M,Nfft,SegNo,N);                  % initialize
Rxx = zeros(M,M,Nfft,N);                    % initialize

accurate_flag = 1;
%% loop over segments
for segIdx = 0:SegNo - 1
    
    T = (segIdx*JUMP + 1):(segIdx*JUMP + Nfft); % current time frame
    x_buf1 = x(:,T,:).*repmat(win,M,1,N); % multiply in analysis window
    x_buf2 = reshape(permute(x_buf1,[2,4,1,3]),[Nfft,1,M,N]);
    x_buf3 = permute(sum(x_buf2,2),[3,1,4,2]);
    X_buf = fft(x_buf3,Nfft,2);
    X(:,:,segIdx + 1,:)= X_buf;
    
end

%% Periodogram power spectral density estimate
for f = 1:Nfft
    normValue = 1/((win*win')*SegNo);
    if accurate_flag
        for nn = 1:N
        Xtemp = squeeze(permute(X(:,f,:,nn),[1,3,2])); % (M x SegNo)
        Rxx(:,:,f,nn) = normValue*(Xtemp*Xtemp');
        end
    else
        Xtemp = squeeze(permute(X(:,f,:,:),[1,3,4,2])); % (M x SegNo)
        Rxx(:,:,f,:) = normValue*(multiprod(Xtemp,permute(conj(Xtemp),[2,1,3])));
    end
end