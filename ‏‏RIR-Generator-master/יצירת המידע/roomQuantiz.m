function [ room_rir_list_quant ] = roomQuantiz( room_rir_list )
% roomQuantiz takes a full room experimnet and quantizes it

partition = 5;

rir_sz = size(room_rir_list);
room_rir_list_quant = cell(rir_sz);

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
            room_rir_list_quant{vlvl,spsp,bb} = cell(size(room_rir_list{vlvl,spsp,bb}));
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
            
            [num_spkrPos,num_of_dis] = size( room_rir_list_quant{vlvl,spsp,bb});
            for ss = 1:num_spkrPos
                fprintf('speaker#: %d/%d \t', ss,num_spkrPos);
                tim = tic;

                for dd = num_of_dis:-1:1
                    if ~isempty(room_rir_list{vlvl,spsp,bb}{ss,dd})
                        if  dd == num_of_dis && ss == 1
                            [room_rir_list_quant{vlvl,spsp,bb}{ss,dd},partition,codebook] = roomRirQunatiz( room_rir_list{vlvl,spsp,bb}{ss,dd},partition);
                        else
                            [room_rir_list_quant{vlvl,spsp,bb}{ss,dd},partition,codebook] = roomRirQunatiz( room_rir_list{vlvl,spsp,bb}{ss,dd},partition,codebook);
                        end
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

