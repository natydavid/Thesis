function [ peak_dens ] = peaks_density( rirs_env )
% peaks_density returns peaks density calculated by number of peaks todivided by time

%% parameters
global strtFromPeak
global Ts

[smp_num,dim]=size(rirs_env);
start_win = 0.1 ; % [s]
win_len = ceil(start_win/Ts);

%% Initialization 
peak_dens = zeros(smp_num,1);

%%
for ss = 1:smp_num
    thisRir = rirs_env(ss,:);
    if strtFromPeak
        [~,peakIdx] = findpeaks(thisRir,'NPeaks',1);
%         [~,peakIdx] = findpeaks(thisRir,'MinPeakWidth',1,'NPeaks',1);
        thisRir = thisRir(peakIdx:end);
        thisRir = thisRir(1:win_len);
    else
        thisRir = thisRir(1:win_len);
    end
    [peaks,~] = findpeaks(thisRir,'MinPeakWidth',1);
    peak_dens(ss) = length(peaks)/start_win;

end
end

