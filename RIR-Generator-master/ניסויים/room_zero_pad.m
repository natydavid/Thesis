function [ room_rirs ] = room_zero_pad( room_rirs )
% room_zero_pad pads all rirs to have the same length

rooms_num = size(room_rirs,1);
room_sz = size(room_rirs{1});
speaker_num = room_sz(1);
CD_num = room_sz(2);

for rr = 1:rooms_num
    dim(rr) = size(room_rirs{rr}{1,1},2);
end

max_dim = max(dim);
pad_room = find(dim~=max_dim);
for rr = pad_room
    for ss = 1:speaker_num
        for cc = 1:CD_num
            sz = size(room_rirs{rr}{ss,cc});
            pad = zeros(sz(1),max_dim-sz(2));
            room_rirs{rr}{ss,cc} = [room_rirs{rr}{ss,cc} pad];
        end
    end
end