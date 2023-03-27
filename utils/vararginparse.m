function vararginparse(cellarray, reqvarstoparse, opvarstoparse)
    if(mod(length(cellarray),2))
        error('Mismatched number of arg names and arg values!')
    end

    if(~isempty(cellarray))
        inargnames = cellarray(1:2:end);
    else
        inargnames = {};
    end

    for i = 1:length(reqvarstoparse)
        if any(strcmp(inargnames, reqvarstoparse{i}))
            ind = find(strcmp(inargnames, reqvarstoparse{i}));
            assignin('caller', cellarray{2*ind-1}, cellarray{2*ind});
        else
            error(sprintf("%s is required but was'nt found", reqvarstoparse{i}))
        end
    end

    for i = 1:length(opvarstoparse)
        if any(strcmp(inargnames, opvarstoparse{i}))
            ind = find(strcmp(inargnames, opvarstoparse{i}));
            assignin('caller', cellarray{2*ind-1}, cellarray{2*ind});
        end
    end
end