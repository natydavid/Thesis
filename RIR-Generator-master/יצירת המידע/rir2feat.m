function [ room_feat_list ] = rir2feat( room_rir_list,save_name,room_list,nworkers )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
rir_sz = size(room_rir_list);
room_feat_list = cell(rir_sz);

vol_cls = rir_sz(1);
shp_cls = rir_sz(2);
if length(rir_sz)==3
    beta_num = rir_sz(3);
else
    beta_num = 1;
end



for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            room_feat_list{vlvl,spsp,bb} = cell(size(room_rir_list{vlvl,spsp,bb}));
        end
    end
end

%tracking parameters
totalIter = vol_cls*shp_cls*beta_num;
n=0;

for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            n=n+1;
            fprintf('Volume#: %d/%d\t Shape#: %d/%d\t Beta#: %d/%d\n',...
                vlvl,vol_cls,...
                spsp,shp_cls,...
                bb,beta_num);
            
            [num_spkrPos,num_of_dis] = size( room_feat_list{vlvl,spsp,bb});
            for ss = 1:num_spkrPos
                fprintf('speaker#: %d/%d \t', ss,num_spkrPos);
                tim = tic;
                
                if exist('nworkers','var')
                    thisRirs  =room_rir_list{vlvl,spsp,bb}(ss,:);
                    temp = cell(1,num_of_dis);
                    parfor (dd = 1:num_of_dis,nworkers)
                        if ~isempty(thisRirs{dd})
                            temp{dd} = rir_feature_extractor( thisRirs{dd});
                        end
                    end
                    room_feat_list{vlvl,spsp,bb}(ss,:) = temp;
                else
                    for dd = 1:num_of_dis
                        if ~isempty(room_rir_list{vlvl,spsp,bb}{ss,dd})
                            room_feat_list{vlvl,spsp,bb}{ss,dd} = rir_feature_extractor( room_rir_list{vlvl,spsp,bb}{ss,dd});
                        end
                    end
                end             
                            
                tim = toc(tim);
                fprintf('[elaps = %.2f s]\n',tim);
            end
            
            prec = (n/totalIter)*100;
            fprintf('Progress: %d%s\n',prec,'%');
        end
    end
end



if exist('save_name','var')
    save(save_name,'room_feat_list');
    room_file = matfile(save_name,'Writable',true);
    room_file.room_list = room_list;
end

end

