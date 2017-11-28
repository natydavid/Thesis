function [ room_rirs ] = room_squeeze( room_rirs )
% room_squeeze squeezes all rirs to have the same length as min

rooms_num = size(room_rirs,1);
room_sz = size(room_rirs{1});
speaker_num = room_sz(1);
CD_num = room_sz(2);

for rr = 1:rooms_num
    dim(rr) = size(room_rirs{rr}{1,1},2);
end

global feat_flag speech_RTF_flag
if feat_flag || speech_RTF_flag
    min_factor = 1;
else
    min_factor = 0.5;
end

if min_factor==1;
    min_dim = min(dim);
    squeeze_room = find(dim~=min_dim);
else
    min_dim = floor(min_factor*min(dim));
    squeeze_room = 1:rooms_num;
end


for rr = squeeze_room
    for ss = 1:speaker_num
        for cc = 1:CD_num
            if ~isempty(room_rirs{rr}{ss,cc})
                room_rirs{rr}{ss,cc} = room_rirs{rr}{ss,cc}(:,1:min_dim);
            end
        end
    end
end