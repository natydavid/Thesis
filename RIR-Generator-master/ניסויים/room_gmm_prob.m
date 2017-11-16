function [ prob_mat ] = room_gmm_prob( room_feat,gmm_models )
% room_gmm_prob recieves a list of gmm models and observation point and retrns
% a matrix (observation X models) with the log porbability of each point to belong to the model
model_num = length(gmm_models);
nPoints = size(room_feat,1);

prob_mat = nan(nPoints,model_num);
for mm = 1:model_num
    prob_mat(:,mm) = gaussmixp(room_feat,gmm_models{mm}.mu,gmm_models{mm}.sigma,gmm_models{mm}.w);
end

end

