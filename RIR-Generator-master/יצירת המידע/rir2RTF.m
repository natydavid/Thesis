function [ room_RTF_list ] = rir2RTF( room_rir_list,room_list,save_name,nworkers )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
rir_sz = size(room_rir_list);
room_RTF_list = cell(rir_sz);

vol_cls = rir_sz(1);
shp_cls = rir_sz(2);
if length(rir_sz)==3
    beta_num = rir_sz(3);
else
    beta_num = 1;
end
fs = 16000;                 % Sample frequency (samples/s)


for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            room_RTF_list{vlvl,spsp,bb} = cell(size(room_rir_list{vlvl,spsp,bb}));
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
            
            [num_spkrPos,num_of_dis] = size( room_RTF_list{vlvl,spsp,bb});
            RT60_len = size( room_rir_list{vlvl,spsp,bb}{1,1},2);
            for ss = 1:num_spkrPos
                fprintf('speaker#: %d/%d \t', ss,num_spkrPos);
                tim = tic;
                
                if exist('nworkers','var')
                    thisRirs  = room_rir_list{vlvl,spsp,bb}(ss,:);
                    temp_rir_est = cell(1,num_of_dis);
                    parfor (dd = 1:num_of_dis,nworkers)
                        if ~isempty(thisRirs{dd})
                            mic1 = thisRirs{dd}(1:2:end,:);
                            mic2 = thisRirs{dd}(2:2:end,:);
                            temp_rir_est{dd} = room_RTF_est( mic1,mic2,fs,RT60_len );
                        end
                    end
                    room_RTF_list{vlvl,spsp,bb}(ss,:) = temp_rir_est;
                else
                    for dd = 1:num_of_dis
                        if ~isempty(room_rir_list{vlvl,spsp,bb}{ss,dd})
                            mic1 = room_rir_list{vlvl,spsp,bb}{ss,dd}(1:2:end,:);
                            mic2 = room_rir_list{vlvl,spsp,bb}{ss,dd}(2:2:end,:);
                            room_RTF_list{vlvl,spsp,bb}{ss,dd} = room_RTF_est( mic1,mic2,fs,RT60_len );
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
    save(save_name,'room_RTF_list');
    room_file = matfile(save_name,'Writable',true);
    room_file.room_list = room_list;
    room_file.fs = fs;
end

end

