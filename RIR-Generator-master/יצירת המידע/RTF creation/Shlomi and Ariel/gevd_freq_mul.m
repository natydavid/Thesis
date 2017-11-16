function [vv,d] = gevd_freq_mul(R,B)

% Generalized Eigenvalue Decomposition (frequency)

% inputs:
% -------
% Rdesd:    signal PSD (M x M x Nfft)
% Rstat:    noise   PSD (M x M x Nfft)

% outputs:
% --------
% vv:       eigen vectors (M x M x Nfft)
% d:        eigen values  (M x Nfft)

% initializations:
% ----------------
M = size(R,1); % number of microphones
Nfft = size(R,3); % fft length
N = size(R,4); % % number of samples
vv = zeros(M,M,Nfft,N); % initialize eigen vectors
d = zeros(M,Nfft,N); % initialize eigen values

% loop over fft bins
% ------------------
for f = 1:Nfft
    Bt = squeeze(B(:,:,f,:));
    Bt = 0.5*(Bt+permute(conj(Bt),[2,1,3]));
    if (f == 1)
        Bt = repmat(eye(M),1,1,N);
    end
    Rt = squeeze(R(:, :, f,:));
    Rt = 0.5*(Rt+permute(conj(Rt),[2,1,3]));
    for nn = 1:N
        sqB = chol(Bt(:,:,nn)); % Cholesky factorization (sqB'*sqB = Bt)
        isqB = inv(sqB); % % inverse sqB
        [Vz,Dz] = eig(isqB'*Rt(:,:,nn)*isqB);
        % Vz: right eigenvectors
        % Dz: Eigenvalues (returned as a diagonal matrix)
        dt = diag(Dz); % Eigenvalues (returned as vector)
        Vz = sqB'*Vz;
        [~,II] = sort(abs(dt),'descend');
        d(:,f,nn) = dt(II);
        vv(:,:,f,nn) = Vz(:,II);
        for e = 1:M
            vv(:,e,f,nn) = vv(:,e,f,nn)/norm(squeeze(vv(:,e,f,nn)));
        end
    end
end

end % function gevd_freq


