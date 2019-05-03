%###############################################################################
% UpdateAMPLData
%
% Determine probabilities of Y sequences and send to AMPL
%###############################################################################
function [] = UpdateAMPLData(ampl, Settings, Data)
    N = size(Data.Y, 1);
    ampl.getParameter('N').setValues(N);

    aYHAT = ampl.getSet('YHAT');
    YHat = cell2mat(cell(aYHAT.get().toArray()));

    YInt = WideToBinary(Data.Y);
    PMF = tabulate(YInt);
    PMF = PMF(:,1:2);

    [C IP ID] = intersect(YHat(:,1), PMF(:,1), 'stable');
    PMF = PMF(ID,1:2);
    assert(length(PMF(:,1)) == length(YHat));
    PMF(:,2) = PMF(:,2)/N;

    vQ = ampl.getParameter('Q');
    Idx = [PMF(:,1)];
    vQ.setValues(Idx, PMF(:,2));
end
