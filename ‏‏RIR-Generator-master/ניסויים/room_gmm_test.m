function room_gmm_test( vargin )
% room_gmm_prob recieves a list of gmm models ,observation point and labels
% and prints a confusion matrix of the results

room_feat_list = vargin.room_rir_list;
room_list  = vargin.room_list;
gmm_models = vargin.gmm_models;
exp_disc = vargin.exp_disc;
speaker = vargin.speaker;
squeeze_pad_flag = vargin.squeeze_pad_flag;
CD = vargin.CD;
log_flag = vargin.log_flag;

% logging
rir_generator_path = which('rir_generator_path.m');
[rir_generator_path,~,~] = fileparts(rir_generator_path);
log_filename = [rir_generator_path '\ניסויים' '\GMM Results.txt'];

room_num = length(gmm_models);

%% squeeze/pad
if squeeze_pad_flag == 1
    room_feat_list = room_squeeze( room_feat_list );
elseif squeeze_pad_flag == 2
    room_feat_list = room_zero_pad( room_feat_list );
end




%% transfering all data to one matrix to deliver to algorithm
room_matrixs = cell(room_num,1);

for rr = 1:room_num
    for ss = speaker
        for cc = CD
            if isfield(vargin,'RTF_peek_samples')
                RTF_peek_samples = vargin.RTF_peek_samples;
                thisRirs = room_feat_list{rr}{ss,cc};
                
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
                
                if isempty(room_matrixs{rr})
                    room_matrixs{rr} = temp_matrix;
                else
                    room_matrixs{rr} = [room_matrixs{rr}; temp_matrix];
                end
            else
                if isempty(room_matrixs{rr})
                    
                    room_matrixs{rr} = room_feat_list{rr}{ss,cc};
                else
                    room_matrixs{rr} = [room_matrixs{rr}; room_feat_list{rr}{ss,cc}];
                end
            end
        end
    end
end


%% calculating log probality for each room and success precantage

conf_mat = nan(room_num,room_num+1);
col_name = cell(1,room_num+1);

for rr = 1:room_num
    if isfield(vargin,'normalize_data')
        if vargin.normalize_data
            min_val = gmm_models{rr}.min_val;
            max_val = gmm_models{rr}.max_val;
            room_matrixs{rr} = normalizeFeature(room_matrixs{rr},min_val,max_val);
        end
    end
    room_prob_matrix = room_gmm_prob( room_matrixs{rr},gmm_models );
    conf_mat(rr,:) = prob2conf(room_prob_matrix,rr);
    col_name{rr} = ['room_' num2str(rr)];
end
roomDisc = extractfield(room_list, 'roomDisc');
col_name{rr+1} = 'success_prec';


Table = myTable(conf_mat,col_name,roomDisc)
if log_flag
    logRoomExpResults( Table,log_filename,exp_disc )
end

end


function [conf_row] = prob2conf(room_prob_matrix,gold_class)
mat_sz = size(room_prob_matrix);
max_room = zeros(mat_sz);
[~,max_idx] = max(room_prob_matrix,[],2);
lin_idx = sub2ind(mat_sz,(1:mat_sz(1))',max_idx);
max_room(lin_idx) = 1;
class_count = sum(max_room,1);
success_prec = (class_count(gold_class)/mat_sz(1))*100;
conf_row = [class_count success_prec];
end
