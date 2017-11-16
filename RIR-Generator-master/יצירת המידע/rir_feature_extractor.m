function [ rir_feature ] = rir_feature_extractor( rirs )

global strtFromPeak
global fs
global Ts
global kurtosis_flag

strtFromPeak = 1;
fs = 16000;
Ts = 1/16000;

kurtosis_flag = 0;

rir_feature = [];

[ Hilbert_env, rms_env  ] = rir_envelope( rirs);

%% Feature extraction - notes uses article notation

% peaks density
rir_feature = [rir_feature peaks_density( rms_env )]; %N_S
rir_feature = [rir_feature peaks_density( Hilbert_env )]; %N_H

% Kurtosis of rir, rms_env, Hilbert_env
rir_feature = [rir_feature kurtosis(rirs,kurtosis_flag,2)]; %K_h
rir_feature = [rir_feature kurtosis(rms_env,kurtosis_flag,2)]; %K_S
rir_feature = [rir_feature kurtosis(Hilbert_env,kurtosis_flag,2)]; %K_S

% frequency features - kurtosis of RTF, std of RTF, number of modes, number of mean crossing
rir_feature = [rir_feature frequency_features( rirs )]; 



end