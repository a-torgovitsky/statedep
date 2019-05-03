%*******************************************************************************
% RunSIPP
%*******************************************************************************
function [] = RunSIPP(SaveDir, SimSet, SimNum, ExitOnEnd)

if ~exist('SaveDir', 'var')
    SaveDir = '';
end

SCRIPTTAG = 'sipp-estimates';
[CleanUpPath CleanUpLog SaveDir] = Setup(SCRIPTTAG, SaveDir);

if ~exist('SimSet', 'var') & ~exist('SimNum', 'var')
    FlagRunWholeCycle = 1;
else
    FlagRunWholeCycle = 0;
end

if ~FlagRunWholeCycle
    [Settings] = LoadSpec(SimSet, SimNum);
    ResultsSubdir = ...
        fullfile(SimSet, sprintf('%03d', SimNum));
    ExecuteThenRecord(Settings, SaveDir, ResultsSubdir);
else
    ThisSimSet = 'main';
    ThisSimNum = 1;

    while ~isempty(ThisSimSet)
        [Settings NextSimSet NextSimNum] = LoadSpec(ThisSimSet, ThisSimNum);
        ResultsSubdir = ...
            fullfile(ThisSimSet, sprintf('%03d', ThisSimNum));

        ExecuteThenRecord(Settings, SaveDir, ResultsSubdir);

        FlagBuildTable = isempty(NextSimSet);
        if FlagBuildTable
            BuildTable(SaveDir, fullfile(SaveDir, 'results', ThisSimSet));
        end
        ThisSimSet = NextSimSet;
        ThisSimNum = NextSimNum;
    end
end

if ~exist('ExitOnEnd', 'var')
    ExitOnEnd = 0;
end
if ExitOnEnd
    exit;
end

end

%*******************************************************************************
%*******************************************************************************
%*******************************************************************************
%*******************************************************************************
% LoadSpec
%*******************************************************************************
%*******************************************************************************
%*******************************************************************************
%*******************************************************************************
function [Settings NextSimSet NextSimNum] = LoadSpec(SimSet, SimNum)

%*******************************************************************************
% General settings
%*******************************************************************************
Settings.T = 6;
Settings.B = 250;
Settings.DataPath = fullfile(pwd, '../data/sipp08.tsv');

Settings.PDBR = 0;
Settings.OptPeriod = (Settings.T+1);
Settings.Tau = .25;
Settings.Parameters = ...
   {'TSD', 'PSD', 'PSD_G0', 'PSD_G00', 'PSD_G1', 'PSD_G11'};
Settings.BuildConfidenceRegions = 1;
Settings.RunMisspecificationTest = 1;
Settings.ParametersToTest = Settings.Parameters;

MaxDimST = Settings.T - 2;
SIGMALIST = 0:.05:.4;

FlagEnd = 0;
switch SimSet
case 'main'
    AfterThisSet = 'sigma';

    switch SimNum

    case 1 % EE only

    case 2
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(0, MaxDimST);

    case 3
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(1, MaxDimST);

    case 4
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(2, MaxDimST);

    case 5
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(3, MaxDimST);

    case 6
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(4, MaxDimST);

    case 7
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(4, MaxDimST);
        Settings.Assumption_MTS = 1;
        Settings.Assumption_DimMTS = 2;

    case 8
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(2, MaxDimST);
        Settings.Assumption_MTR = 1;

    case 9
        Settings.Assumption_TIV = 1;

    % PDBR below
    case 10
        Settings.PDBR = 1;
        Settings.TimeDummies = 0;

    case 11
        Settings.PDBR = 1;
        Settings.TimeDummies = 1;

    case 12
        Settings.PDBR = 1;
        Settings.TimeDummies = 1;
        Settings.Age = 1;
        FlagEnd = 1;

    otherwise
        error('SimNum not recognized.')

    end

case 'sigma'
    AfterThisSet = 'sigma-young';

    Settings.Parameters = ...
       {'PSD_G0', 'PSD_G00'};
    Settings.ParametersToTest = Settings.Parameters;
    Settings.CalculateMaxImpliedChange = 1;

    Settings.Assumption_ST = 1;
    Settings.Assumption_DimST = MaxDimST;
    if SimNum <= length(SIGMALIST)
        Settings.Assumption_SigmaST = SIGMALIST(SimNum);
    else
        error('SimNum not recognized.')
    end
    FlagEnd = (SimNum == length(SIGMALIST));

case 'sigma-young'
    AfterThisSet = '';
    Settings.DataPath = '../data/sipp08-young.tsv';

    Settings.Parameters = ...
       {'PSD_G0', 'PSD_G00'};
    Settings.ParametersToTest = Settings.Parameters;

    Settings.Assumption_ST = 1;
    Settings.Assumption_DimST = MaxDimST;
    if SimNum <= length(SIGMALIST)
        Settings.Assumption_SigmaST = SIGMALIST(SimNum);
    else
        error('SimNum not recognized.')
    end
    FlagEnd = (SimNum == length(SIGMALIST));

case 'extra'
    AfterThisSet = '';

    switch SimNum

    case 1 % Empirical evidence only

    case 2 % MATR only
        Settings.Assumption_MATR = 1;

    case 3 % ST(2) only
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(2, MaxDimST);

    case 4 % ST(2)+ DSC
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(2, MaxDimST);
        Settings.Assumption_DSC = 1;

    case 5 % ST(2) + MATR
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(2, MaxDimST);
        Settings.Assumption_MATR = 1;

    case 6 % ST(2) + DSC + MATR
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(2, MaxDimST);
        Settings.Assumption_DSC = 1;
        Settings.Assumption_MATR = 1;

    case 7 % ST(2) + MTS(2)
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(2, MaxDimST);
        Settings.Assumption_MTS = 1;
        Settings.Assumption_DimMTS = 2;

    case 8 % ST(2) + MTS(2) + DSC + MATR
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = min(2, MaxDimST);
        Settings.Assumption_MTS = 1;
        Settings.Assumption_DimMTS = 2;
        Settings.Assumption_DSC = 1;
        Settings.Assumption_MATR = 1;

    otherwise
        error('SimNum not recognized.')

    end
otherwise
    error('SimSet not recognized.')
end

if FlagEnd
    NextSimNum = 1;
    NextSimSet = AfterThisSet;
else
    NextSimSet = SimSet;
    NextSimNum = SimNum + 1;
end

end

%*******************************************************************************
%*******************************************************************************
%*******************************************************************************
%*******************************************************************************
%*******************************************************************************
%*******************************************************************************
%*******************************************************************************
function ExecuteThenRecord( Settings,...
                            SaveDir,...
                            ResultsSubdir)
%*******************************************************************************
%*******************************************************************************
    OriginalPath = CreateResultsDir(SaveDir, ResultsSubdir);

    RecordStructure(Settings, 'SettingsBefore.out');
    if ~Settings.PDBR
        Settings = rmfield(Settings, 'PDBR');
        [Results Settings] = DPO(Settings);
    else
        Settings = rmfield(Settings, 'PDBR');
        fid = fopen('PDBR.out', 'w'); % helps w/ table script
        fprintf(fid, '1');
        fclose(fid);
        [Results Settings] = PDBR(Settings);
    end

    RecordStructure(Settings, 'SettingsAfter.out');
    RecordAssumptions(Results.AssumptionString, 'Assumptions.out');
    RecordBounds(Settings.Parameters, Results.Bounds, 'Bounds.out')

    RecordSingleNumber( Results.MinCriterion,...
                        'MinCriterion',...
                        'MinCriterion.out');

    if isfield(Settings, 'CalculateMaxImpliedChange')
        dlmwrite('MaxImpliedChange.out', Results.MaxImpliedChange,...
                 'delimiter', '\t',...
                 'precision', 5);
    end

    if Settings.BuildConfidenceRegions
        RecordCR(Settings, Results, 'ConfidenceRegions');
    end

    RecordSingleNumber( Results.MSPValue,...
                        'p-value',...
                        'Misspecification.out');

    if isfield(Results, 'CFHN')
        RecordStructure(Results.CFHN, 'CFHN.out');
    end

    cd(OriginalPath);
end

function RecordCR(Settings, Results, OutfilenameStub)
    for j = 1:1:length(Settings.Tests)
        for a = 1:1:length(Settings.LevelsCR)
            fid = fopen([OutfilenameStub '_'...
                'A' int2str(Settings.LevelsCR(a)*100)...
                '_' Settings.Tests{j} '.out'], 'w');
            for p = 1:1:length(Settings.ParametersToTest)
                fprintf(fid, '%-15s %12.10f %12.10f\n',...
                    Settings.ParametersToTest{p},...
                    Results.CR(p,j,a,1), Results.CR(p,j,a,2));
            end
        end
    end
end

function [] = BuildTable(SaveDir, ResultsDir)
    cd(fullfile(SaveDir, 'post'));
    pwd
    cmd = ['python BuildTableBounds.py ' ResultsDir]
    system(cmd);
end
