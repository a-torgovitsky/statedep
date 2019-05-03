%*******************************************************************************
% ResampleData
%
% Draw a random sample of size S with or without replacement from Data.
% Data is a structure with each field corresponding to a variable.
%*******************************************************************************
function DataResample = ResampleData(Data, S, replacement, Seed)
    if ~(S >= 0)
        error('Resample size is not a positive integer.');
    end
    if ~(Seed >= 0)
        error('Seed is not a positive integer.');
    end

    rng(Seed); % The right way
    %rand('seed', Seed); % The wrong way

    FieldNames = fieldnames(Data);
    N = size(Data.(FieldNames{1}), 1);
    I = randsample(N, S, replacement);

    for f = 1:1:length(FieldNames)
        DIM = ndims(Data.(FieldNames{f}));
        if (DIM == 2)
            DataResample.(FieldNames{f}) = Data.(FieldNames{f})(I,:);
        elseif (DIM == 3)
            DataResample.(FieldNames{f}) = Data.(FieldNames{f})(I,:,:);
        else
            error('DIM is %d -- that''s strange!', DIM);
        end
    end
end
