function dna = rna2dna2(rna)
% rna2dna2 converts an RNA sequence into a DNA sequence.
%Unlike the standard Matlab function, this ignores if the sequence has
%non-GCAUT bases

if ~ischar(rna)
    dna = rna;
    return
end

dna = strrep(rna,'U','T');
dna = strrep(dna,'u','t');
