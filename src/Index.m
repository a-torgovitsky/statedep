%*******************************************************************************
% Index
%
% List is a cell array of strings, e.g. {'Test1', 'Test2', 'Test5'}
% So Index('Test5', List) would give me back 3.
%*******************************************************************************
function i = Index(Item, List)
    [tf i] = ismember(Item, List);

    if ~tf
        error('Item not in List');
    end
end
