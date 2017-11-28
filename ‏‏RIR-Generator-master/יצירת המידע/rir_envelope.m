function [ Hilbert_env, rms_env  ] = rir_envelope( rirs)
% rir_envelope returns the upper and lower envelope of the rir using 
% the Hillbert transform method and rms method

%% paramters

[smp_num,dim]=size(rirs);
rms_win = 3;
filter_len = 4; % even number returns a smother envelope

%% Initalizations 

Hilbert_env = zeros(smp_num,dim);
rms_env = zeros(smp_num,dim);



%%
for ss = 1:smp_num
    
    [rms_env(ss,:),~] = envelope(rirs(ss,:),rms_win,'rms');    
    [Hilbert_env(ss,:),~] = envelope(rirs(ss,:),filter_len,'analytic');
    
end
end

