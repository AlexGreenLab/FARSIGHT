function i = issame2(x,y)
% function i = issame2(x,y)

i = 0;
if isequal(double(x),double(y))
    if ~max(isnan(x)) && ~max(isnan(y))
        i = 1;
        if ischar(x) == ischar(y)
            i = 2;
        end
    end
end
