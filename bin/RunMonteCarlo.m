%*******************************************************************************
% RunMonteCarlo
%*******************************************************************************
function RunMonteCarlo(SaveDir, SimNumber, NMultiplier)
if ~exist('SaveDir', 'var')
    SaveDir = '';
end
if ~exist('SimNumber')
    error('Must provide SimNumber.')
end
%*******************************************************************************
%*******************************************************************************
%*******************************************************************************
SCRIPTTAG = 'mc';
[CleanUpPath CleanUpLog SaveDir] = Setup(SCRIPTTAG, SaveDir);
%*******************************************************************************
%*******************************************************************************
%*******************************************************************************
Settings.DataPath = fullfile(pwd, '../data/sipp08.tsv');

MCSettings.M = 500;

switch SimNumber
    case 1 % Point estimates T = 6, ST(4)
        Settings.T = 6;
        Settings.OptPeriod = Settings.T + 1;
        Settings.Parameters = {'PSD', 'PSD_G0', 'PSD_G00', 'PSD_G1', 'PSD_G11'};

        Settings.TestAListOfPoints = 0;

        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = Settings.T - 2;

    case 2 % Confidence regions
        Settings.T = 3;
        Settings.OptPeriod = Settings.T + 1;
        Settings.Parameters = {'PSD_G0'};

        Settings.TestAListOfPoints = 1;
        Settings.B = 500;
        Settings.Tests = {'CNS', 'SS'};
        Settings.PointsToTest{1} = [.09 .15 .21 .27 .30 .48 .51 .54 .60 .66];

        Settings.Assumption_MTR = 1;
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = 1;

    case 3 % Point estimates, T = 3, ST(1) and MTR -- like the CR one
        Settings.T = 3;
        Settings.OptPeriod = Settings.T + 1;
        Settings.Parameters = {'PSD', 'PSD_G0', 'PSD_G00', 'PSD_G1', 'PSD_G11'};

        Settings.TestAListOfPoints = 0;

        Settings.Assumption_MTR = 1;
        Settings.Assumption_ST = 1;
        Settings.Assumption_DimST = 1;

    otherwise
        error('SimNumber not recognized.');
end

if ~exist('NMultiplier')
    MCSettings.NMultiplier = 1;
    DirName = fullfile(SaveDir, 'results');
else
    MCSettings.NMultiplier = NMultiplier;
    DirName = fullfile( SaveDir, 'results',...
                        sprintf('N%.1f', MCSettings.NMultiplier));
    mkdir(DirName);
end

cd(DirName);
MonteCarlo(Settings, MCSettings);
