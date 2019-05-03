%*******************************************************************************
% EstimateIdentifiedSet
%
% Estimate the minimum value of the criterion (or optionally assume it is 0).
% Then estimate the identified set by max'ing and min'ing the target parameter
% subject to the criterion being smaller than this value (x (1 + Tau)).
% Return those bounds and the minimum value of the criterion.
%
% Bounds will be +/- Inf if infeasible, but this is only possible if the sample
% identified set is empty and the user opts to assume that it is non-empty.
%*******************************************************************************
function [Bounds MinCriterion] = EstimateIdentifiedSet(ampl, Settings)
    if (Settings.Noise >= 1)
        disp('Finding minimum criterion ...');
    end

    if Settings.AssumeIDSetNonempty
        MinCriterion = 0;

        if (Settings.Noise >= 1)
            disp(sprintf('\tAssumed to be 0.'));
        end
    else
        % This basically re-uses the misspecification test code
        %   what you put for ``Points'' in ComputeTestStatistics doesn't matter
        %   as long as it is scalar
        Settings.ActiveParam = 'MS';
        MinCriterion = ComputeTestStatistics(ampl, [-123], Settings);

        if (Settings.Noise >= 1)
            str = sprintf('\tComputed to be %19.16f.', MinCriterion);
            disp(str);
            disp(repmat('=', 1, Settings.DisplaySepLen));
        end
    end

    Bounds = zeros(length(Settings.Parameters), 2);
    Bounds(:,1) = +Inf;
    Bounds(:,2) = -Inf;

    if MinCriterion < Settings.DeclareCriterionToBeZeroTol
        if (Settings.Noise >= 1)
            disp('Computing sample identified sets...')
        end
        ChangeOptimizationProblem(ampl, Settings, 'IdentifiedSet');
        FlagCompute = 1;
    else
        if (Settings.Noise >= 1)
            disp('Estimating identified sets...')
        end

        ChangeOptimizationProblem(ampl, Settings, 'Estimation');
        ampl.getParameter('CriterionHat').setValues(MinCriterion);
        FlagCompute = 0;
    end

    for j = 1:1:length(Settings.Parameters)
        ObjectiveName = ['min' Settings.Parameters{j}];
        ampl.eval(['objective ' ObjectiveName ';']);


        if Settings.NoisyOptimization
            eval('ampl.solve');
        else
            evalc('ampl.solve');
        end
        ErrorCheckOptimization(ampl,...
            ['EstimateIdentifiedSet: ', ObjectiveName], 0);
        SolveResult = ampl.getValue('solve_result');

        if (strcmp(SolveResult, 'infeasible'))
            FlagInfeasible = 1;
        else
            FlagInfeasible = 0;
        end

        if FlagCompute & FlagInfeasible % Return all bounds as +Inf/-Inf
            return;
        end

        if ~FlagCompute & FlagInfeasible
            error(['Estimation problem is infeasible.'...
                   ' This should not happen.'])
        end

        % If we passed both of the above checks, record and proceed to ub
        Bounds(j,1) = SafelyGetObjective(ampl, ObjectiveName);
        ObjectiveName = ['max' Settings.Parameters{j}];
        ampl.eval(['objective ' ObjectiveName ';']);
        evalc('ampl.solve');
        ErrorCheckOptimization(ampl,...
            ['EstimateIdentifiedSet: ' ObjectiveName], 0);
        SolveResult = ampl.getValue('solve_result');
        assert(strcmp(SolveResult, 'infeasible') == 0);
        Bounds(j,2) = SafelyGetObjective(ampl, ObjectiveName);
    end
end
