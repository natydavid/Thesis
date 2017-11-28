function [ features ] = frequency_features( rirs )
% frequency_features extract features related to frequency:
% k_rtf - kurtosis of RTF
% s_rtf - std of RTF
% n_rtf - number of modes
% mc_rtf - number of mean crossing

[smp_num,dim]=size(rirs);

global kurtosis_flag
global fs

f0 = 100; %[Hz]
k{1} = 1:3;
k{2} = 1:5;
k{3} = 1:10;
dev_factor = [2 1.5 1.2];
dev_num = length(dev_factor);

feat_num = 4;
band_num = length([k{1} k{2} k{3}]);
features_tmp = zeros(smp_num,band_num,feat_num);

n = 2^nextpow2(dim);
freq_axe = ((1:n)-1)*(fs/n);

rtf = abs(fft(rirs,n,2));

FC = 0; % band counter
for dd = 1:dev_num
    for kk = k{dd}
        
        FC = FC+1;
        
        freq_start = (k{dd}(kk)*f0)/dev_factor(dd);
        freq_end = k{dd}(kk)*f0*dev_factor(dd);
        seg_start = find(freq_axe<=freq_start,1,'last');
        seg_end = find(freq_axe>=freq_end,1);
        rtf_seg = rtf(:,seg_start:seg_end);
        
        k_rtf = kurtosis(rtf_seg,kurtosis_flag,2);
        s_rtf = std(rtf_seg,0,2);
        n_rtf = findpeaks_mat(rtf_seg,smp_num);
        mc_rtf = mean_crossing(rtf_seg);
        
        features_tmp(:,FC,1) = k_rtf;
        features_tmp(:,FC,2) = s_rtf;
        features_tmp(:,FC,3) = n_rtf;
        features_tmp(:,FC,4) = mc_rtf;
        
        
    end
end

features = [];
for ff = 1:feat_num
    features = [features features_tmp(:,:,ff)];
end


end

function [mc_rtf] = mean_crossing(rtf_seg)
mu = mean(rtf_seg,2);
mc_logic = bsxfun(@gt, rtf_seg, mu);
mc_rtf = sum(mc_logic,2);

end

function [n_rtf] = findpeaks_mat(rtf_seg,smp_num)
n_rtf = zeros(smp_num,1);
for ss = 1:smp_num
    thisRtf = rtf_seg(ss,:);
    [peaks,~] = findpeaks(thisRtf,'MinPeakWidth',1);
    n_rtf(ss) = length(peaks);
    
end
end
