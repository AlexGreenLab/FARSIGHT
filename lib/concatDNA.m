function [sequence,bases_removed] = concatDNA(varargin)
%function [sequence,bases_removed] = concatDNA(varargin)

DEBUG = 0;
sequence = varargin{1};
bases_removed = 0;

for c1 = 1:length(varargin)
    varargin{c1} = rna2dna2(varargin{c1});
end

for c1 = 2:length(varargin)
    start_base = 1;
    prev = varargin{c1-1};
    curr = varargin{c1};
    start_index = 1;
    for start_base = 1:min(length(prev),length(curr))
        if issame2(prev(end-start_base+1:end), curr(1:start_base)) == 2
            start_index = start_base + 1;
        end
    end
    bases_removed = bases_removed + start_index - 1;

    sequence = [sequence,curr(start_index:end)];
end

if DEBUG == 1
    acc = 0;
    for c1 = 1:length(varargin)
        temp = length(findstr(sequence,varargin{c1}));
        acc = temp + acc;
        if temp > 1
            fprintf('\t%d duplicates for part %s.\n',temp-1,varargin{c1});
        end
    end
    fprintf('Acc = %d, length(varargin) = %d\n',acc,length(varargin));
end