function [ room_feat_list ] = rir2feat( room_rir_list,save_name )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
rir_sz = size(room_rir_list);
room_feat_list = cell(rir_sz);

vol_cls = rir_sz(1);
shp_cls = rir_sz(2);
beta_num = rir_sz(3);



for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            room_feat_list{vlvl,spsp,bb} = cell(size(room_feat_list{vlvl,spsp,bb}));
        end
    end
end


for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            
            [num_spkrPos,num_of_dis] = size( room_feat_list{vlvl,spsp,bb});
            for ss = 1:num_spkrPos
                for dd = 1:num_of_dis
                    room_feat_list{vlvl,spsp,bb}{ss,dd} = rir_feature_extractor( room_feat_list{vlvl,spsp,bb}{ss,dd});
                end
            end

        end
    end
end

if exist('save_name','var')
    save_name = [save_name '.mat'];
    save(save_name,'room_feat_list');
end

end

