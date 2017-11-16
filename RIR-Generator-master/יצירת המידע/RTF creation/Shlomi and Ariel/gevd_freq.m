function [vv,d] = gevd_freq(R,B)

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
vv = zeros(M,M,Nfft); % initialize eigen vectors
d = zeros(M,Nfft); % initialize eigen values

% loop over fft bins
% ------------------
for f = 1:Nfft
    Bt = squeeze(B(:,:,f));
    Bt = 0.5*(Bt+Bt');
    if (f == 1)
        Bt = eye(M);
    end
    Rt = squeeze(R(:, :, f));
    Rt = 0.5*(Rt+Rt');
    sqB = chol(Bt); % Cholesky factorization (sqB'*sqB = Bt)
    isqB = inv(sqB); % % inverse sqB
    [Vz,Dz] = eig(isqB'*Rt*isqB);
    % Vz: right eigenvectors
    % Dz: Eigenvalues (returned as a diagonal matrix)
    dt = diag(Dz); % Eigenvalues (returned as vector)
    Vz = sqB'*Vz;
    [~,II] = sort(abs(dt),'descend');
    d(:,f) = dt(II);
    vv(:,:,f) = Vz(:,II);   
    for e = 1:M
        vv(:,e,f) = vv(:,e,f)/norm(squeeze(vv(:,e,f)));
    end
end

end % function gevd_freq


