%###############################################################################
% ComputeTestStatistic
%
% Compute the sample test statistic at every point in Points.
%###############################################################################
function [TS] = ComputeTestStatistics(ampl, Points, Settings)

    TS = zeros(length(Points), 1);

    ChangeOptimizationProblem(ampl, Settings, 'Criterion');

    % Solve at each point in Points
    ampl.getConstraint('FixParameter').restore;
    for t = 1:1:length(Points)
        ampl.getParameter('Fix').setValues(Points(t));
        if Settings.NoisyOptimization
            eval('ampl.solve');
        else
            evalc('ampl.solve');
        end

        IDStr = ['ComputeTestStatistics: '];
        if isfield(Settings, 'MCPoints')
            IDStr = [IDStr 'm = ' int2str(Settings.CurrentSim) ', '];
        end
        IDStr = [IDStr ' t = ' num2str(Points(t))];
        ErrorCheckOptimization(ampl, IDStr, 1);

        TS(t,1) = SafelyGetObjective(ampl, 'minCriterion');
    end
end
