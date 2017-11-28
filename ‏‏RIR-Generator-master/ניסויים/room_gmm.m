function [ gmm_models ] = room_gmm( vargin )
% this function receives rooms' rirs and creates a gmm model for each room
% room_rirs is a rooms X 1 cell array
% room is speakerPointsInRoom X CD cell array
% input:
% varin.room_rir_list -  room rirs
% varin.nmix - number of mixtures for gmm
% varin.speaker - a vector indicating which speaker points should be used for training
% varin.CD - a vector indicating which CD to take
% varin.save_name - name to save the model
%
% output:
% varout.gmm_models -  rooms X 1 cell array each cell contains a gmm model

room_rirs = vargin.room_rir_list;
nmix = vargin.nmix;
speaker = vargin.speaker;
CD = vargin.CD ;
squeeze_pad_flag = vargin.squeeze_pad_flag;
min_samples = inf;

rir_sz = size(room_rirs);

vol_cls = rir_sz(1);
shp_cls = rir_sz(2);
if length(rir_sz)==3
    beta_num = rir_sz(3);
else
    beta_num = 1;
end

room_sz = size(room_rirs{1});
speaker_num = room_sz(1);
CD_num = room_sz(2);

%% swithcing dim and number of samples to fit gmm_em function
if 0
    
    for rr = 1:rooms_num
        for ss = 1:speaker_num
            for cc = 1:CD_num
                room_rirs{rr}{ss,cc} = room_rirs{rr}{ss,cc}';
            end
        end
    end
end


%% squeeze/pad
if squeeze_pad_flag == 1
    room_rirs = room_squeeze( room_rirs );
elseif squeeze_pad_flag == 2
    room_rirs = room_zero_pad( room_rirs );
end


%% transfering all data to one matrix to deliver to algorithm
room_matrixs = cell(vol_cls,shp_cls,beta_num);

for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            
            for ss = speaker
                for cc = CD
                    if isfield(vargin,'RTF_peek_samples')
                        RTF_peek_samples = vargin.RTF_peek_samples;
                        thisRirs = room_rirs{vlvl,spsp,bb}{ss,cc};
                        
                        [~,idx] = max(thisRirs,[],2);
                        rowSub = repmat((1:size(thisRirs,1))',1,RTF_peek_samples*2+1);
                        colSub = ones(size(thisRirs,1),RTF_peek_samples*2+1);
                        colSub(:,1) = idx-RTF_peek_samples;
                        colSub = cumsum(colSub,2);
                        takenMsk = false(size(thisRirs));
                        mskIdx = sub2ind(size(thisRirs),rowSub(:),colSub(:));
                        takenMsk(mskIdx) = true;
                        thisRirs_tr = thisRirs';
                        takenMsk_tr = takenMsk';
                        temp = thisRirs_tr(takenMsk_tr);
                        temp = reshape(temp,RTF_peek_samples*2+1,size(thisRirs,1));
                        temp_matrix = temp';
                        
                        if isempty(room_matrixs{vlvl,spsp,bb})
                            room_matrixs{vlvl,spsp,bb} = temp_matrix;
                        else
                            room_matrixs{vlvl,spsp,bb} = [room_matrixs{vlvl,spsp,bb}; temp_matrix];
                        end
                    else
                        if isempty(room_matrixs{vlvl,spsp,bb})
                            room_matrixs{vlvl,spsp,bb} = room_rirs{vlvl,spsp,bb}{ss,cc};
                        else
                            room_matrixs{vlvl,spsp,bb} = [room_matrixs{vlvl,spsp,bb}; room_rirs{vlvl,spsp,bb}{ss,cc}];
                        end
                    end
                end
            end
            min_samples = min(min_samples,size(room_matrixs{vlvl,spsp,bb},1));
        end
    end
end


if 0
    
    %% reducing to the same number of samples for each class
    for vlvl = 1:vol_cls
        for spsp = 1:shp_cls
            for bb = 1:beta_num
                num_of_samples = size(room_matrixs{vlvl,spsp,bb},1);
                if num_of_samples > min_samples
                    idx_taken = unidrnd(num_of_samples,[1 min_samples]);
                    room_matrixs{vlvl,spsp,bb} = room_matrixs{vlvl,spsp,bb}(idx_taken,:);
                end
            end
        end
    end
end

%% normalize data
if isfield(vargin,'normalize_data')
    if vargin.normalize_data
        total_data = cell2mat(room_matrixs(:));
        vargout = normalizeFeature( total_data );
        min_val = vargout.min_val;
        max_val = vargout.max_val;
        
        for vlvl = 1:vol_cls
            for spsp = 1:shp_cls
                for bb = 1:beta_num
                    room_matrixs{vlvl,spsp,bb} = normalizeFeature(room_matrixs{vlvl,spsp,bb},min_val,max_val);
                end
            end
        end
        
    end
end

%% training models

gmm = struct('w',nan,...
    'mu',nan,...
    'sigma',nan...
    );

if isfield(vargin,'normalize_data')
    if vargin.normalize_data
        gmm.min_val = min_val;
        gmm.max_val = max_val;
    end
end

gmm_models = repmat({gmm},vol_cls,shp_cls,beta_num);

for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            [gmm.mu, gmm.sigma, gmm.w] = gaussmix(room_matrixs{vlvl,spsp,bb},[],[],nmix);
            gmm_models{vlvl,spsp,bb} = gmm ;
        end
    end
end

%% saving models
if isfield(vargin,'save_name')
    save_name = vargin.save_name ;
    save(save_name,'gmm_models');
end
end

