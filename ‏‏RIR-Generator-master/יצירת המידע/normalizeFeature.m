function [ vargout ] = normalizeFeature( feature,min_val,max_val )
[nSmp,~] = size(feature);

if exist('min_val','var') %for case you already know the norm parameters
    min_rep = repmat(min_val,nSmp,1);
    feature = feature - min_rep;
    max_val(max_val==0) = 1;
    max_rep = repmat(max_val,nSmp,1);
    feature = feature./max_rep;
    vargout = feature;
else
    min_val = min(feature,[],1);
    min_rep = repmat(min_val,nSmp,1);
    feature = feature - min_rep;
    max_val = max(feature,[],1);
    max_val(max_val==0) = 1;
    max_rep = repmat(max_val,nSmp,1);
    feature = feature./max_rep;
    vargout.feature = feature;
    vargout.min_val = min_val;
    vargout.max_val = max_val;
end


end

