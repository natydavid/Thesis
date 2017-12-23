function [ diffuse_noise ] = diffusionNoise( C, n, nworkers )
% Inputs:
% 
% C - STFT mixing matrix
% n - M signals in the time domain [L x M x nSamples]

% Outputs:
% 
% C - STFT mixing matrix

nSamples = size(n,3); % Number of sanmples
M = size(n,2); % Number of sensors
L = size(n,1); % Length input signal
K = (size(C,3)-1)*2;
temp = stft(n(:,1,1), K, K/4, 1).';
win_size = size(temp,2);

% Initialization
N = zeros(M,win_size,K,nSamples);
X = zeros(size(N));  % STFT output matrix
X(:,:,1,:) = X(1,1,1,1);

% Short-time Fourier transform

for m = 1 : M
    tempSignal = squeeze(n(:,m,:));
    parfor (nn = 1:nSamples,nworkers)
        N(m,:,:,nn) = stft(tempSignal(:,nn), K, K/4, 1).';
    end
end

% Generate output in the STFT domain for each frequency bin k
for k = 2:K/2+1
    tempN = squeeze(N(:,:,k,:));
    temp = multiprod(C(:,:,k)' * tempN);
    X(:,:,k,:) = permute(temp,[1 2 4 3]);
end

% Inverse STFT
temp = real(istft(squeeze(X(m,:,:,1)).', K, K/4, 1));
diffuse_noise = zeros(size(temp,1),M,nSamples);
for m = 1 : M
    tempSignal = squeeze(X(m,:,:,:));
    parfor (nn = 1:nSamples,nworkers)
        diffuse_noise(:,m,nn) = real(istft(tempSignal(:,:,nn).', K, K/4, 1));
    end
end
diffuse_noise = diffuse_noise(1:L,:,:)

end



