%###############################################################################
% Create a string that summarizes the assumptions, and output to the user
%   The string gets returned through the Results struct as well for recording
%###############################################################################
function S = CreateAssumptionString(Settings)

    Names = {'MTR'};
    Values = {sprintf('%d', Settings.Assumption_MTR)};

    Names = [Names ['MATR']];
    Values = [Values sprintf('%d', Settings.Assumption_MATR)];

    Names = [Names ['ST']];
    Values = [Values sprintf('%d', Settings.Assumption_ST)];

    Names = [Names ['mST']];
    Values = [Values sprintf('%d',...
        Settings.Assumption_ST*Settings.Assumption_DimST...
        + ...
        (1-Settings.Assumption_ST)*-1)];

    Names = [Names ['SigmaST']];
    Values = [Values sprintf('%f',...
        Settings.Assumption_ST*Settings.Assumption_SigmaST...
        + ...
        (1-Settings.Assumption_ST)*-1)];

    Names = [Names ['TIV']];
    Values = [Values sprintf('%d', Settings.Assumption_TIV)];

    Names = [Names ['DSC']];
    Values = [Values sprintf('%d', Settings.Assumption_DSC)];

    Names = [Names ['MIV']];
    Values = [Values sprintf('%d', Settings.Assumption_MIV)];

    Names = [Names ['MTS']];
    Values = [Values sprintf('%d', Settings.Assumption_MTS)];

    Names = [Names ['qMTS']];
    Values = [Values sprintf('%d', ...
        Settings.Assumption_MTS*Settings.Assumption_DimMTS...
        + ...
        (1-Settings.Assumption_MTS)*-1)];

    format = '%-15s %-s';
    for j = 1:1:length(Names)
        S{j} = sprintf(format, Names{j}, Values{j});
    end
end
