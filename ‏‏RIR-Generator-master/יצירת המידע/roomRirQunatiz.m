function [ room_rir_quant,partition,codebook ] = roomRirQunatiz( room_rir,partition,codebook )
% roomRirQunatiz quntize a matrix using the lloyd algorithm to first find a good partition and a codebook

method = 'round';% 'round' 'quant' 'simple'

[nSamples, dim] = size(room_rir);


switch method
    
    case 'quant'
        room_rir = room_rir'; % in order to work correct in reshape
        room_rir_reshp = reshape(room_rir,[],1);
        if nargin<3 % first time - guessing a partition and a code book
            N = partition;
            min_val = min(room_rir_reshp);
            max_val = max(room_rir_reshp);
            gap = (max_val-min_val)/(N-1);
            codebook = min_val:gap:max_val;
            [partition,codebook] = lloyds(room_rir_reshp,codebook);
        end
        
        [~,room_rir_quant_reshp] = quantiz(room_rir_reshp,partition,codebook);
        room_rir_quant_reshp = room_rir_quant_reshp';
        room_rir_quant = reshape(room_rir_quant_reshp,dim,nSamples)';
        
    case 'simple'
        
        if floor(partition) == partition
            N = partition;
            min_val = min(min(room_rir));
            max_val = max(max(room_rir));
            res = (max_val-min_val)/(N-1);
            partition = res;
        else
            res = partition;
        end
        room_rir_quant = quant(room_rir,res);
        codebook = [];
        
    case 'round'
        
        dec_point = partition;
        room_rir_quant = round(room_rir,dec_point);
        codebook = [];
end

end

