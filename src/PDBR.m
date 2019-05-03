%###############################################################################
% PDBR
%
% Estimate a parametric dynamic binary response (PDBR) model
% In particular, a dynamic random effects probit
%###############################################################################
function [Results Settings] = PDBR(SettingsIn)
%###############################################################################
% Define defaults, apply input
%###############################################################################
Settings.DataPath = '';
Settings.T = [];
Settings.TimeDummies = 1;
Settings.Age = 0;
Settings.Parameters = { 'ATE', 'PSD', 'NSD',...
                        'PSD_G0', 'PSD_G00',...
                        'PSD_G1', 'PSD_G11'};

Settings.InitialSeed = 2;
Settings.B = 250;
Settings.LevelsCR = .05; % (oops -- this is actually 1 - level)

Settings.NumGHPoints = 16;
Settings.Noise = 0;
Settings.BSProgressFrequency = 20;
Settings.OptParamsFile = 'opt-params.csv';
Settings.TgtParamsFile = 'tgt-params.csv';

if ~isstruct(SettingsIn)
    error('Expected structure for input. Quitting.');
else
    Settings = UpdateStruct(Settings, SettingsIn, 0, 0);
end

%###############################################################################
% LOAD DATA
%###############################################################################
[Settings Data] = LoadData(Settings);

%###############################################################################
% ADD TIME DUMMIES
%
% Note -- the other covariate is initial period outcomes
% But this not treated as part of "X" the way I coded it (see PDBR.mod)
%###############################################################################
N = size(Data.Y, 1);
if Settings.TimeDummies
    % Time dummies
    Data.X = zeros(N, Settings.T, Settings.T);
    for t = 1:1:Settings.T
        Data.X(:,t,t) = ones(N,1);
    end
else % just include a constant term in this case
    Data.X = ones(N, Settings.T, 1);
end

if Settings.Age
    K = size(Data.X, 3);
    for t = 1:1:Settings.T
        Data.X(:,t, K + 1) = Data.Age(:,t);
    end
end

%###############################################################################
% CREATE AMPL INSTANCE
%###############################################################################
ampl = InitializeAMPL({'PDBR.mod'}, Settings);
CleanUpAMPL = onCleanup(@()ampl.close());
ampl.setOption('solver', 'knitro');

% Quadrature stuff for AMPL -- only needs to be done once
[nodes weights] = GaussHermite(Settings.NumGHPoints);
ampl.getParameter('J').setValues(Settings.NumGHPoints);
ampl.getParameter('ghnodes').setValues((1:1:Settings.NumGHPoints)', nodes);
ampl.getParameter('ghweights').setValues((1:1:Settings.NumGHPoints)', weights);

%###############################################################################
% POINT ESTIMATES
%###############################################################################
SendDataToAMPL(ampl, Data);
MaximizeLikelihood(ampl, Settings.Noise);
OptT = UpdateOptParamTable(ampl);
TgtT = UpdateTgtParamTable(ampl, Settings.Parameters);

Results.TgtParam = TgtT;
Results.OptParam = OptT;

Results.Bounds = repmat(Results.TgtParam{1,:}', 1, length(Settings.Parameters));
Results.MinCriterion = 0;
Results.MSPValue = -1;

%###############################################################################
% BOOTSTRAP
%###############################################################################
if Settings.B <= 0
    return;
else
    Settings.BuildConfidenceRegions = 1;
end

for b = 1:1:Settings.B
    if (mod(b, Settings.BSProgressFrequency) == 0)
        disp(sprintf('Starting bootstrap replication %d.', b));
        TgtRes
        OptRes
    end
    DataBS = ResampleData(Data, N, 1, Settings.InitialSeed + b);
    SendDataToAMPL(ampl, DataBS)
    MaximizeLikelihood(ampl, Settings.Noise);
    OptT = UpdateOptParamTable(ampl, OptT);
    TgtT = UpdateTgtParamTable(ampl, Settings.Parameters, TgtT);

    [TgtRes, OptRes] = ConfidenceIntervals(TgtT, OptT, Settings.LevelsCR(1));

    writetable(OptRes, Settings.OptParamsFile, 'WriteRowNames', true);
    writetable(TgtRes, Settings.TgtParamsFile, 'WriteRowNames', true);
end

Results.TgtParam = TgtRes;
Results.OptParam = OptRes;

% For consistency with the DPO code
Results.CR = zeros(length(Settings.Parameters), 1, 1, 2);
Results.CR(:,1,1,1) = TgtRes{1,:}';
Results.CR(:,1,1,2) = TgtRes{3,:}';
Settings.Tests = {'CNS'}; % meaningless -- just determines file name
Settings.ParametersToTest = Settings.Parameters;
Settings.Assumption_MTR = 0;
Settings.Assumption_MATR = 0;
if Settings.TimeDummies
    Settings.Assumption_ST = 0;
    Settings.Assumption_DimST = 0;
else
    Settings.Assumption_ST = 1;
    Settings.Assumption_DimST = 4;
end
Settings.Assumption_SigmaST = 0;
Settings.Assumption_TIV = 1;
Settings.Assumption_DSC = 0; % is actually true, but lets leave it out
Settings.Assumption_MIV = 0;
Settings.Assumption_MTS = 1;
Settings.Assumption_DimMTS = 6;
Results.AssumptionString = CreateAssumptionString(Settings);

end

%###############################################################################
%###############################################################################
%###############################################################################
%###############################################################################
%###############################################################################
%###############################################################################
%###############################################################################
%###############################################################################
%###############################################################################
%###############################################################################
%###############################################################################
%###############################################################################
function SendDataToAMPL(ampl, Data)

    N = size(Data.Y, 1);
    T = size(Data.Y, 2) - 1;
    K = size(Data.X, 3);

    ampl.getParameter('N').setValues(N);
    ampl.getParameter('T').setValues(T);
    ampl.getParameter('K').setValues(K);

    tic;

    Idx = [];
    Val = [];
    for t = 0:1:T
        Idx = [Idx; (1:1:N)' t*ones(N,1)];
        Val = [Val; Data.Y(:,t+1)];
    end
    ampl.getParameter('y').setValues(Idx, Val);

    Idx = [];
    Val = [];
    for t = 1:1:T
        for k = 1:1:K
            Idx = [Idx; (1:1:N)' t*ones(N,1) k*ones(N,1)];
            Val = [Val; Data.X(:,t,k)];
        end
    end
    ampl.getParameter('x').setValues(Idx, Val);

    disp(sprintf('Initializing data took %5.3f seconds', toc));
end

%###############################################################################
%###############################################################################
%###############################################################################

function SolveTime = MaximizeLikelihood(ampl, Noise)
    if ~exist('Noise')
        Noise = 0;
    end

    tic;
    if Noise
        eval('ampl.solve');
    else
        evalc('ampl.solve');
    end
    SolveTime = toc;

    if Noise
        disp(sprintf('Solving took %5.3f seconds.', toc));
    end

    SolveResult = ampl.getValue('solve_result');
    if ~ismember(SolveResult, {'solved', 'solved?'})
        error('Sample likelihood solve_result was %s', SolveResult);
    end
end

%###############################################################################
%###############################################################################
%###############################################################################

function OptParamTable = UpdateOptParamTable(ampl, OptParamTable)
    z = ampl.getObjective('LLH');
    llh = z.value;

    gamma = ampl.getValue('gamma');
    lambda = ampl.getValue('lambda');
    sigma = ampl.getValue('sigma');

    AddTable = table(llh, gamma, lambda, sigma);

    v = ampl.getVariable('beta');
    vDF = v.getValues('val');
    beta = vDF.getColumnAsDoubles('val');

    for k = 1:1:length(beta)
        AddTable.(sprintf('beta%d', k)) = beta(k);
    end

    if ~exist('OptParamTable')
        OptParamTable = AddTable;
    else
        OptParamTable = [OptParamTable; AddTable];
    end
end

%###############################################################################
%###############################################################################
%###############################################################################

function TgtParamTable = UpdateTgtParamTable(ampl, ParamList, TgtParamTable)
    params = nan(1, length(ParamList));
    for p = 1:1:length(ParamList)
        params(p) = ampl.getValue(strcat(ParamList{p}, '[T+1]'));
    end
    AddTable = array2table(params);
    AddTable.Properties.VariableNames = ParamList;

    if ~exist('TgtParamTable')
        TgtParamTable = AddTable;
    else
        TgtParamTable = [TgtParamTable; AddTable];
    end
end

%###############################################################################
%###############################################################################
%###############################################################################
function [TgtRes, OptRes] = ConfidenceIntervals(TgtT, OptT, Level)

    OptRes = repmat(OptT(1,:), 3, 1);
    lowq = Level/2;
    upq = 1 - lowq;
    for q = 1:1:width(OptRes)
        OptRes{1,q} = ComputeQuantile(table2array(OptT(2:end,q)), lowq);
        OptRes{3,q} = ComputeQuantile(table2array(OptT(2:end,q)), upq);
    end
    OptRes.Properties.RowNames = {'lb', 'est', 'ub'};

    TgtRes = repmat(TgtT(1,:), 3, 1);
    for q = 1:1:width(TgtRes)
        TgtRes{1,q} = ComputeQuantile(table2array(TgtT(2:end,q)), lowq);
        TgtRes{3,q} = ComputeQuantile(table2array(TgtT(2:end,q)), upq);
    end
    TgtRes.Properties.RowNames = {'lb', 'est', 'ub'};
end

%###############################################################################
%###############################################################################
%###############################################################################
%###############################################################################
% 04/18/19
% https://groups.google.com/d/msg/ampl/mZBnmtvXQSE/UqudPWhQCwAJ
%###############################################################################
function t = DfToTable(df)
    t = table;
    h = df.getHeaders();
    for i=1:length(h)
        %t.(strrep(char(h(i)), '.', '_')) = df.getColumn(h(i));
        t.(strrep(char(h(i)), '.', '_')) = cell2mat(df.getColumn(h(i)));
    end
end
