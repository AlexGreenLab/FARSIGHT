function rna = dna2rna2(dna)
%dna2rna2 converts a DNA sequence into an RNA sequence.
%Unlike the standard Matlab function, this ignores if the sequence has
%non-GCAUT bases

if ~ischar(dna)
    rna = dna;
    return
end

rna = strrep(dna,'T','U');
rna = strrep(rna,'t','u');