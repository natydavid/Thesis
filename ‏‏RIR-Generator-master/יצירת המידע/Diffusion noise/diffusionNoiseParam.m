function [ C ] = diffusionNoiseParam( type_nf, Fs, K, M, d, method )
% Inputs:
% 
% type_nf -  Type of noise field: 'spherical' or 'cylindrical'
% Fs - Sample frequency (Hz)
% K - FFT length
% M - Number of sensors
% d - Inter sensor distance (m)
% method - 'cholesky' or 'eigen'
% 
% Outputs:
% C - STFT mixing matrix

%% Initialization
c = 340;                  % Sound velocity (m/s)
C = zeros(M,M,K/2+1); % STFT mixing matrix


%% Generate matrix with desired spatial coherence
ww = 2*pi*Fs*(0:K/2)/K;
DC = zeros(M,M,K/2+1);
for p = 1:M
    for q = 1:M
        if p == q
            DC(p,q,:) = ones(1,1,K/2+1);
        else
            switch lower(type_nf)
                case 'spherical'
                    DC(p,q,:) = sinc(ww*abs(p-q)*d/(c*pi));
                    
                case 'cylindrical'
                    DC(p,q,:) = bessel(0,ww*abs(p-q)*d/c);
                    
                otherwise
                    error('Unknown noise field.')
            end
        end
    end
end

% Generate output in the STFT domain for each frequency bin k
for k = 2:K/2+1
    switch lower(method)
        case 'cholesky'
            C(:,:,k) = chol(DC(:,:,k));
            
        case 'eigen'
            [V,D] = eig(DC(:,:,k));
            C(:,:,k) = sqrt(D) * V';
            
        otherwise
            error('Unknown method specified.');
    end
    
end

end

