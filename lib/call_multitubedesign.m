function result = call_multitubedesign(design_file_name,prevent_array,check_indices,IS_DNA)
% function result = call_multitubedesign(design_file_name,prevent_array,check_indices,IS_DNA)

[status,result] = system(sprintf('multitubedesign %s',design_file_name));
if isempty(result)
    if IS_DNA
        [prevent_string0,~] = generatePreventDNASeq(prevent_array,'A');
    else
        [prevent_string0,~] = generatePreventRNASeq(prevent_array,'A');
    end
    seq_output = parseNUPACKfile2SeqInfo(sprintf('%s_0.npo',design_file_name(1:end-3)));
    for c1 = 1:length(check_indices)
        if IS_DNA
            [temp_prevent_string,~] = generatePreventDNASeq(prevent_array,seq_output{check_indices(c1),2});
        else
            [temp_prevent_string,~] = generatePreventRNASeq(prevent_array,seq_output{check_indices(c1),2});
        end
        %temp_prevent_string
        %prevent_array
        %prevent_string0
        %fprintf('ERROR: Prevent requirements not met for %s.\n',seq_output{check_indices(c1),2});
        if ~issame2(temp_prevent_string,prevent_string0)
            result = 'ERROR: Prevent requirements not met.';
            return;
        end
    end
end


% [prevent_string,prevent_array] = generatePreventDNASeq(base_prevented,target_seq)
% 
% temp_prevent_array = [];
% for cx2 = 1:size(prevent_array0,1)
%     temp_prevent_array(cx2,1:2) = sprintf('%s%d',prevent_array0{cx2}(1),length(prevent_array0{cx2})+prevent_extra);
% end
% temp_prevent_array = char(temp_prevent_array);
% 
% 
