function output_string = d_revcomp(input_string)
%function output_string = d_revcomp(input_string)
%Any non-ACGUT bases will be ignored and replaced with Ns

input_string = rna2dna2(input_string);
output_string = char('N'+zeros(1,length(input_string)));
for c1 = 1:length(input_string)
    if input_string(c1) == 'A'
        output_string(c1) = 'T';
    elseif input_string(c1) == 'C'
        output_string(c1) = 'G';
    elseif input_string(c1) == 'G'
        output_string(c1) = 'C';
    elseif input_string(c1) == 'T'
        output_string(c1) = 'A';
    end
end
output_string = fliplr(output_string);
