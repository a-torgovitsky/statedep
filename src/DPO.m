%###############################################################################
% DPO
%
% Construct identified sets and/or confidence regions for the DPO model.
% Many parameters may be specified, but only one set of assumptions can
% be used for each call.
%
% Input:
%   SettingsIn
%       a structure containing all options (see defaults below)
%
% Output:
%   Results
%       a structure of results
%   Settings
%       the updated structure of settings, which might be useful
%   Data
%       the data -- useful if running a Monte Carlo
%###############################################################################
function [Results Settings Data] = DPO(SettingsIn, DataIn)
%###############################################################################
% Define defaults, apply input
%###############################################################################
% Settings that control the data, model, modelling assumptions, and parameters
Settings.DataPath = '';
Settings.T = [];
Settings.OptPeriod = [];
Settings.Parameters = {'PSD', 'PSD_G0', 'PSD_G00', 'PSD_G1', 'PSD_G11'};
Settings.Assumption_MTR = 0;
Settings.Assumption_MATR = 0;
Settings.Assumption_ST = 0;
Settings.Assumption_DimST = 0;
Settings.Assumption_SigmaST = 0;
Settings.Assumption_TIV = 0;
Settings.Assumption_DSC = 0;
Settings.Assumption_MIV = 0;
Settings.Assumption_MTS = 0;
Settings.Assumption_DimMTS = 2;

% Settings that control what procedures are done
Settings.BuildConfidenceRegions = 0;
Settings.RunMisspecificationTest = 0;
Settings.TestAListOfPoints = 0;
Settings.ParametersToTest = {};
Settings.PointsToTest = {};
Settings.AssumeIDSetNonempty = 0;
Settings.ComputeCFHNBounds = 1;
Settings.CalculateMaxImpliedChange = 1;

% Settings that control options or tuning parameters for statistical inference
Settings.B = 500;
Settings.InitialSeed = 1;
Settings.Tau = .25;
Settings.SSExp = 2/3;
Settings.Tests = {'CNS'};
Settings.LevelsCR = [.05];
Settings.LevelsTestList = [.01 .05 .10];

% Less important numerical tuning parameters and solver options
Settings.Solver = 'cplex';
Settings.BracketTol = 1e-3;
Settings.SkipTestingTol = 1e-7;
Settings.RejectTol = 1e-6;
Settings.FeasTolDefault = 1e-6; % CPLEX default
Settings.FeasTolStepFactor = 10;
Settings.FeasTolMax = 1e-3;
Settings.PreSolveEps = 1e-10;
Settings.DeclareCriterionToBeZeroTol = 1e-6;

% Output options
Settings.Noise = 1;
Settings.NoisyOptimization = 0;
Settings.DisplaySepLen = 80;
Settings.GetDefaultSettings = 0; % Just replace default settings then return

if ~isstruct(SettingsIn)
    error('Expected structure for input. Quitting.');
else
    Settings = UpdateStruct(Settings, SettingsIn, 1);
end

if Settings.GetDefaultSettings
    Results = [];
    Data = [];
    Settings.GetDefaultSettings = 0;
    return;
end

%###############################################################################
% Some error checking on parameter specification
%###############################################################################
ParameterUniverse = ...
    {'ATE', 'PSD', ...
     'PSD_G0', 'PSD_G00', ...
     'PSD_G1', 'PSD_G11',...
     'NSD', 'TSD'};
if ~all(ismember(Settings.Parameters, ParameterUniverse))
    error('Settings.Parameters is incorrectly specified.')
end
if (Settings.TestAListOfPoints | Settings.BuildConfidenceRegions)
    if ~all(ismember(Settings.ParametersToTest, Settings.Parameters))
        error('Settings.ParametersToTest should be within Settings.Parameters.')
    end
end
if (Settings.TestAListOfPoints)
    if ~(length(Settings.ParametersToTest) == length(Settings.PointsToTest))
        error(['Settings.ParametersToTest and Settings.PointsToTest'...
                ' should be structs of the same length.']);
    end
end

%###############################################################################
% Load data into Matlab (not AMPL)
% If DataIn was passed then bypass this -- this is only used for Monte Carlos
%###############################################################################
if ~exist('DataIn')
    [Settings Data] = LoadData(Settings);
else
    Data = DataIn;
    Settings.N = size(Data.Y, 1);
    assert(Settings.T == (size(Data.Y, 2) - 1));
end

%###############################################################################
% Initialize an instance of AMPL
% Set some solver options (inside InitializeAMPL)
% Create set definitions (this can take a while)
%###############################################################################
ampl = InitializeAMPL({'DPO.mod'}, Settings);
CleanUpAMPL = onCleanup(@()ampl.close());

% Create set definitions in AMPL
CreateAMPLSets(ampl, Settings, Data);

%###############################################################################
% Update data (probabilities) in AMPL
%###############################################################################
UpdateAMPLData(ampl, Settings, Data);
ampl.eval('reset data Q_Sample;')

%###############################################################################
% Summarize assumptions for recording/outputing later
%###############################################################################
Results.AssumptionString = CreateAssumptionString(Settings);
% Output to user
if (Settings.Noise >= 1)
    disp('Beginning computations with the following assumptions:');
    for j = 1:1:length(Results.AssumptionString)
        disp(Results.AssumptionString{j})
    end
    disp(repmat('=', 1, Settings.DisplaySepLen));
end

%###############################################################################
% Calculate max implied change in PSD under Sigma_ST
%###############################################################################
if Settings.CalculateMaxImpliedChange
    Results.MaxImpliedChange = zeros(Settings.T, Settings.T);
    ChangeOptimizationProblem(ampl, Settings, 'MaxImpliedChange');

    for t = 1:1:Settings.T
        for tt = 1:1:Settings.T
            if (tt == t)
                continue; % obviously 0 in this case
            end
            ObjectiveName = sprintf('ChangeInPSD[%d, %d];', t, tt);
            ampl.eval(['objective ' ObjectiveName ';']);

            if Settings.NoisyOptimization
                eval('ampl.solve');
            else
                evalc('ampl.solve');
            end

            SolveResult = ampl.getValue('solve_result');
            assert(~strcmp(SolveResult, 'infeasible')); % should not be infeas

            Results.MaxImpliedChange(t,tt) = ...
                SafelyGetObjective(ampl, ObjectiveName);
        end
    end
else
    Results.MaxImpliedChange = -1*ones(Settings.T, Settings.T);
end

%###############################################################################
% Estimate bounds -- always attempt to do this first
%###############################################################################
[Results.Bounds Results.MinCriterion] = EstimateIdentifiedSet(ampl, Settings);

if (Settings.Noise >= 1)
    DisplayTable = table(transpose(cellstr(Settings.Parameters)),...
                         Results.Bounds(:,1),...
                         Results.Bounds(:,2));

    DisplayTable.Properties.VariableNames = ...
        {'Parameter', 'LB', 'UB'};

    disp('Done. Bounds:')
    disp(DisplayTable);
    disp(repmat('=', 1, Settings.DisplaySepLen));
end

%###############################################################################
% Test a list of points
%###############################################################################
if Settings.TestAListOfPoints
    if ~isnumeric(Settings.LevelsTestList) | ~all(Settings.LevelsTestList > 0)
        error('LevelsTestList is incorrectly specified.')
    end

    for p = 1:1:length(Settings.ParametersToTest)
        Settings.ActiveParam = Settings.ParametersToTest{p};
        if (Settings.Noise >= 1)
            str = sprintf('Testing these points for %s:\n\t',...
                Settings.ActiveParam);
            str = [str sprintf('%5.3f  ', Settings.PointsToTest{p})];
            disp(str);
        end

        [Results.TS{p} Results.PValue{p} Results.CV{p} Results.Reject{p}] ...
            = TestListOfPoints(ampl, Settings, Data,...
                Settings.PointsToTest{p}, Settings.LevelsTestList);

        if (Settings.Noise >= 1)
            DisplayTable = table(Settings.PointsToTest{p}, Results.PValue{p})
            DisplayTable.Properties.VariableNames = {'Point', 'PValue'};

            disp('Done. P-values:')
            disp(DisplayTable);
            disp(repmat('=', 1, Settings.DisplaySepLen));
        end
    end
    if (Settings.Noise >= 1)
        disp(repmat('=', 1, Settings.DisplaySepLen));
    end
else
    Results.TS{1} = -1*ones(1, 1);
    Results.PValue{1} = -1*ones(1, length(Settings.Tests));
    Results.CV{1} = -1*ones(length(Settings.LevelsTestList),...
                            1,...
                            length(Settings.Tests));
    Results.Reject{1} = Results.CV{1};
end

%###############################################################################
% Construct confidence intervals -- if flag is on
%###############################################################################
if Settings.BuildConfidenceRegions
    if ~isnumeric(Settings.LevelsCR) | ~all(Settings.LevelsCR > 0)
        error('LevelsCR is incorrectly specified.')
    end
    if (Settings.Noise >= 1)
        disp('Building confidence regions...')
    end
    Results.CR = BuildConfidenceRegions(ampl, Settings, Data, Results.Bounds);

    if (Settings.Noise >= 1)
        disp(repmat('=', 1, Settings.DisplaySepLen));
    end
else
    % Return dummy info
    Results.CR = nan(  1,...
                       length(Settings.Tests),...
                       length(Settings.LevelsCR),...
                       2);
    CR(:,:,:,1) = 1;
    CR(:,:,:,2) = -1;
end

%###############################################################################
% Conduct misspecification test
%###############################################################################
if (Settings.RunMisspecificationTest)
    if (Settings.Noise >= 1)
        disp('Conducting a misspecification test...')
    end
    Settings.ActiveParam = 'MS';
    % Note that what you put for "Points" is not important w/ MS
    % as long as it is a scalar (otherwise you run it multiple times)
    [Results.MSTS Results.MSPValue] = ...
        TestListOfPoints(ampl, Settings, Data, [-123]);
    if (Settings.Noise >= 1)
        disp(sprintf('p-value was %7.5f', Results.MSPValue));
        disp(repmat('=', 1, Settings.DisplaySepLen));
    end
else
    Results.MSTS = -1;
    Results.MSPValue = -1*ones( 1, length(Settings.Tests));
end

if (Settings.Noise >= 1)
    disp(repmat('=', 1, Settings.DisplaySepLen));
    disp(repmat('=', 1, Settings.DisplaySepLen));
    disp('All done with DPO.')
    disp(repmat('=', 1, Settings.DisplaySepLen));
    disp(repmat('=', 1, Settings.DisplaySepLen));
    disp(repmat('=', 1, Settings.DisplaySepLen));
end

%###############################################################################
% Compute CFHN bounds
%###############################################################################
if Settings.ComputeCFHNBounds
    Results.CFHN = ComputeCFHNBounds(Settings, Data);
else
    Results.CFHN = struct;
end

end
