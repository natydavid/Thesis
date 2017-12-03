function  [data]= python_orgnize_data_no_split( vargin )
% this function receives rooms' rirs and creates a gmm model for each room
% room_rirs is a rooms X 1 cell array
% room is speakerPointsInRoom X CD cell array
% input:
% vargin.room_rir_list -  room rirs
% vargin.train_speaker_pos - a vector indicating which speaker points should be used for training
% vargin.CD - a vector indicating which CD to take
% vargin.save_name - name to save the model
% vargin.save_path - path to save the model


room_rir_list = vargin.room_rir_list;
CD = vargin.CD ;
room_list = vargin.room_list ;
normalize_data = vargin.normalize_data;
if isfield(vargin,'num_of_chunks')
    num_of_chunks = vargin.num_of_chunks;
end
rooms_num = size(room_rir_list,1);
room_sz = size(room_rir_list{1});
speaker_num = room_sz(1);
CD_num = room_sz(2);


%% zero padding the data

room_rir_list = room_zero_pad( room_rir_list );

%% transfering all data to one matrix to deliver to algorithm

room_matrixs = cell(rooms_num,1);
critical_distance_label_rooms = cell(rooms_num,1);
points_num = zeros(1,rooms_num);
for rr = 1:rooms_num
    for ss = 1:speaker_num
        for cc = CD
            if isfield(vargin,'RTF_peek_samples')
                RTF_peek_samples = vargin.RTF_peek_samples;
                thisRirs = room_rir_list{rr}{ss,cc};
                
                [~,idx] = max(thisRirs,[],2);
                rowSub = repmat((1:size(thisRirs,1))',1,RTF_peek_samples*2+1);
                colSub = ones(size(thisRirs,1),RTF_peek_samples*2+1);
                colSub(:,1) = idx-RTF_peek_samples;
                colSub = cumsum(colSub,2);
                takenMsk = false(size(thisRirs));
                
                % throw over range samples
                rowSub = rowSub(:);
                colSub = colSub(:);
                
                taken_idx = colSub <= size(thisRirs,2) & colSub >= 1 ;
                rows_thrown = unique(rowSub(~taken_idx));
                taken_idx = ~ismember(rowSub,rows_thrown); 
                rowSub = rowSub(taken_idx);
                colSub = colSub(taken_idx);
                
                mskIdx = sub2ind(size(thisRirs),rowSub(:),colSub(:));
                takenMsk(mskIdx) = true;
                thisRirs_tr = thisRirs';
                takenMsk_tr = takenMsk';
                temp = thisRirs_tr(takenMsk_tr);
                temp = reshape(temp,RTF_peek_samples*2+1,size(thisRirs,1) - length(rows_thrown));
                temp_matrix = temp';
                
                if isempty(room_matrixs{rr})
                    room_matrixs{rr} = temp_matrix;
                else
                    room_matrixs{rr} = [room_matrixs{rr}; temp_matrix];
                end
                nSamples = size(temp_matrix,1);
            else
                if isempty(room_matrixs{rr})
                    room_matrixs{rr} = room_rir_list{rr}{ss,cc};
                else
                    room_matrixs{rr} = [room_matrixs{rr}; room_rir_list{rr}{ss,cc}];
                end
                nSamples = size(room_rir_list{rr}{ss,cc},1);
            end
            critical_distance_label_rooms{rr} = [critical_distance_label_rooms{rr} ; cc*ones(nSamples,1)];
        end
    end
    points_num(rr) = size(room_matrixs{rr},1);
end




label = [];
feature = [];
critical_distance_label = [];
target_names = {};
if floor(max(vargin.relative_part*points_num))<= min(points_num)
    min_points_num = floor(max(vargin.relative_part*points_num));
else
    min_points_num = floor(vargin.relative_part*min(points_num));
end

for rr = 1:rooms_num
    if isfield(vargin,'same_amount_flag')
        if vargin.same_amount_flag
            taken_idx = unidrnd(points_num(rr),[1 min_points_num]);
            points_num(rr) = min_points_num;
            room_matrixs{rr} = room_matrixs{rr}(taken_idx,:);
            critical_distance_label_rooms{rr} = critical_distance_label_rooms{rr}(taken_idx,:);
        else
            taken_idx = unidrnd(points_num(rr),[1 vargin.relative_part*points_num(rr)]);
            points_num(rr) = vargin.relative_part*points_num(rr);
            room_matrixs{rr} = room_matrixs{rr}(taken_idx,:);
            critical_distance_label_rooms{rr} = critical_distance_label_rooms{rr}(taken_idx,:);
        end
    end
    critical_distance_label = [critical_distance_label ; critical_distance_label_rooms{rr} ];
    label = [label ; (rr-1)*ones(points_num(rr),1)];
    feature = [feature ; room_matrixs{rr}];
    target_names = [target_names room_list(rr).roomDisc];
    
end

if isfield(vargin,'normalize_data')
    if normalize_data
        vargout = normalizeFeature( feature );
        feature = vargout.feature;
        min_val = vargout.min_val;
        max_val = vargout.max_val;
    end
end

points = length(label);

if isfield(vargin,'binary')
    if vargin.binary
        precntile = 20;
        precThresh = prctile(abs(feature),precntile,2);
        %         precThresh = prctile(abs(feature(:)),precntile);
        bin_feature = feature;
        lesOrEq = bsxfun(@le,abs(bin_feature),precThresh);
        bin_feature(lesOrEq)=0;
        grtr = bsxfun(@gt,abs(bin_feature),precThresh);
        bin_feature(grtr)=1;
        data.feature = bin_feature;
    else
        
        data.feature = feature;
    end
else
    data.feature = feature;
end
% data.feature = feature(1:floor(0.1*points),:);
% data.label = label(1:floor(0.1*points));
data.label = label;
data.critical_distance_label = critical_distance_label;
data.points = points;
data.target_names = target_names;


%% saving models

if isfield(vargin,'save_path')
    cd(vargin.save_path);
end

if isfield(vargin,'save_name')
    if num_of_chunks==1
        save_name = [vargin.save_name '.mat'];
        save(save_name,'data');
    end
end


if isfield(vargin,'num_of_chunks')
    if num_of_chunks~=1
        points_per_split = floor(points/num_of_chunks);
        for ss = 1:num_of_chunks
            start_idx = 1+(ss-1)*points_per_split;
            end_idx = ss*points_per_split;
            data.feature = feature(start_idx:end_idx,:);
            data.label = label(start_idx:end_idx);
            data.critical_distance_label = critical_distance_label(start_idx:end_idx);
            data.points = points_per_split;
            save_name = [vargin.save_name '_chunk_' num2str(ss) '.mat'];
            save(save_name,'data');
        end
    end
end
end

