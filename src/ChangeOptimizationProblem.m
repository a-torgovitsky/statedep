%*******************************************************************************
% ChangeOptimizationProblem
%
% Drop/restore assumptions depending on the type of problem to be solved.
%*******************************************************************************
function [] = ChangeOptimizationProblem(ampl, Settings, Type)
    AcceptedTypes = {'IdentifiedSet', 'Estimation', 'Criterion',...
                     'CNS', 'MaxImpliedChange'};
    assert(ismember(Type, AcceptedTypes));

    %***************************************************************************
    % Assumption settings that are applicable for all problem types
    %***************************************************************************
    ampl.getParameter('OptPeriod').setValues(Settings.OptPeriod);

    if ~isnumeric(Settings.Tau) | (Settings.Tau < 0)
        error('Tau is not a positive real number.')
    end
    ampl.getParameter('Tau').setValues(Settings.Tau);

    if Settings.Assumption_MTR
        ampl.getConstraint('MTR').restore;
    else
        ampl.getConstraint('MTR').drop;
    end

    if Settings.Assumption_MATR
        ampl.getConstraint('MATR').restore;
    else
        ampl.getConstraint('MATR').drop;
    end

    if Settings.Assumption_ST
        if (Settings.Assumption_SigmaST == 0)
            ampl.getConstraint('ST').restore;
            ampl.getConstraint('STL').drop;
            ampl.getConstraint('STR').drop;
        else
            ampl.getConstraint('ST').drop;
            ampl.getConstraint('STL').restore;
            ampl.getConstraint('STR').restore;
        end
    else
        ampl.getConstraint('ST').drop;
        ampl.getConstraint('STL').drop;
        ampl.getConstraint('STR').drop;
    end

    if Settings.Assumption_DSC
        ampl.getConstraint('DSC_Forw').restore;
        ampl.getConstraint('DSC_Back').restore;
    else
        ampl.getConstraint('DSC_Forw').drop;
        ampl.getConstraint('DSC_Back').drop;
    end

    if Settings.Assumption_TIV
        ampl.getConstraint('TIV').restore;
    else
        ampl.getConstraint('TIV').drop;
    end

    %***************************************************************************
    % Toggle settings for specific problem types
    %***************************************************************************
    % ObsEq is the exact observational equivalence condition
    %   so this is only imposed when computing the identified set
    if strcmp(Type, 'IdentifiedSet')
        ampl.getConstraint('ObsEq').restore;
    else
        ampl.getConstraint('ObsEq').drop;
    end

    % Fixing parameters is for testing only
    if ismember(Type, {'IdentifiedSet', 'Estimation', 'MaxImpliedChange'})
        ampl.getConstraint('FixParameter').drop;
    else
        ampl.getConstraint('FixParameter').restore;
        if isempty(Settings.ActiveParam)
            error(['Must have Settings.ActiveParam non-empty if Type '...
                   'is not IdentifiedSet.'])
        end
        ampl.getParameter('ActiveParam').setValues(Settings.ActiveParam);
    end

    % MTS is part of the criterion unless you are directly computing
    if ismember(Type, {'IdentifiedSet', 'MaxImpliedChange'})
        if (Settings.Assumption_MTS)
            ampl.getConstraint('MTS').restore;
        else
            ampl.getConstraint('MTS').drop;
        end
    else
        ampl.getConstraint('MTS').drop;
        if (Settings.Assumption_MTS)
            ampl.getParameter('Flag_MTS_In_Objective').setValues(1);
        else
            ampl.getParameter('Flag_MTS_In_Objective').setValues(0);
        end
    end

    % Constrain obs-eq values to lie within a tolerance of the minimum
    % Also require any ``conditional'' parameters to make sense (be in [0,1])
    if strcmp(Type, 'Estimation')
        ampl.getConstraint('EstimatedIDSet').restore;
        ampl.getConstraint('PSD_G0_Bound').restore;
        ampl.getConstraint('PSD_G00_Bound').restore;
        ampl.getConstraint('PSD_G1_Bound').restore;
        ampl.getConstraint('PSD_G11_Bound').restore;
    else
        ampl.getConstraint('EstimatedIDSet').drop;
        ampl.getConstraint('PSD_G0_Bound').drop;
        ampl.getConstraint('PSD_G00_Bound').drop;
        ampl.getConstraint('PSD_G1_Bound').drop;
        ampl.getConstraint('PSD_G11_Bound').drop;
    end

    % A ``_Sample'' version of the above is used in CNS since the probabilities
    % are fixed at their values from the data.
    % There's probably a clean way to combine this with the above, but I thought
    % it was easier like this.
    if strcmp(Type, 'CNS')
        ampl.getConstraint('EstimatedIDSet_Sample').restore;
        ampl.getConstraint('ObsEq_Gap_Sample_Abs').restore;
    else
        ampl.getConstraint('EstimatedIDSet_Sample').drop;
        ampl.getConstraint('ObsEq_Gap_Sample_Abs').drop;
    end

    if strcmp(Type, 'Criterion')
        ampl.eval('objective minCriterion;');
    end

    % These slack variables are used in the criterion function, so they are
    % not needed for directly computing the identified set.
    % For CNS they are replaced by a "Sample" and "Bootstrap" version, so they
    % aren't used there either.
    if ismember(Type, {'Estimation', 'Criterion'})
        ampl.getConstraint('ObsEq_Gap_Abs_Pos').restore;
        ampl.getConstraint('ObsEq_Gap_Abs_Neg').restore;
        ampl.getConstraint('MTS_Gap_Eq_Abs_Pos').restore;
        ampl.getConstraint('MTS_Gap_Eq_Abs_Neg').restore;
    else
        ampl.getConstraint('ObsEq_Gap_Abs_Pos').drop;
        ampl.getConstraint('ObsEq_Gap_Abs_Neg').drop;
        ampl.getConstraint('MTS_Gap_Eq_Abs_Pos').drop;
        ampl.getConstraint('MTS_Gap_Eq_Abs_Neg').drop;
    end

    if strcmp(Type, 'CNS')
        SetTolerance(ampl, Settings.FeasTolDefault);
    end

    % Parameters that are clearly only used in CNS
    if strcmp(Type, 'CNS')
        ampl.eval('objective minCriterion_CNS;');
        ampl.getConstraint('FixParameter_CNS').restore;
        ampl.getConstraint('ZeroOne_CNS').restore;
        ampl.getConstraint('Probability_CNS').restore;

        ampl.getConstraint('ObsEq_Gap_CNS_Abs_Pos').restore;
        ampl.getConstraint('ObsEq_Gap_CNS_Abs_Neg').restore;

        if Settings.Assumption_MTR
            ampl.getConstraint('MTR_CNS').restore;
        else
            ampl.getConstraint('MTR_CNS').drop;
        end

        if Settings.Assumption_MATR
            error('Have not coded MATR for CNS yet');
            %ampl.getConstraint('MATR_CNS').restore;
        end
        %else
            %ampl.getConstraint('MATR_CNS').drop;
        %end

        if Settings.Assumption_ST & (Settings.Assumption_SigmaST == 0)
            ampl.getConstraint('ST_CNS').restore;
        else
            ampl.getConstraint('ST_CNS').drop;
        end

        if Settings.Assumption_ST & (Settings.Assumption_SigmaST > 0)
            ampl.getConstraint('STL_CNS').restore;
            ampl.getConstraint('STR_CNS').restore;
        else
            ampl.getConstraint('STL_CNS').drop;
            ampl.getConstraint('STR_CNS').drop;
        end

        if Settings.Assumption_TIV
            ampl.getConstraint('TIV_CNS').restore;
        else
            ampl.getConstraint('TIV_CNS').drop;
        end

        if Settings.Assumption_MTS
            ampl.getConstraint('MTS_CNS').restore;
            ampl.getConstraint('MTS_Gap_Eq_CNS_Abs_Pos').restore;
            ampl.getConstraint('MTS_Gap_Eq_CNS_Abs_Neg').restore;
            ampl.getConstraint('MTS_Gap_Eq_Sample_Abs').restore;
        else
            ampl.getConstraint('MTS_CNS').drop;
            ampl.getConstraint('MTS_Gap_Eq_CNS_Abs_Pos').drop;
            ampl.getConstraint('MTS_Gap_Eq_CNS_Abs_Neg').drop;
            ampl.getConstraint('MTS_Gap_Eq_Sample_Abs').drop;
        end
    else
        ampl.getConstraint('FixParameter_CNS').drop;
        ampl.getConstraint('ZeroOne_CNS').drop;
        ampl.getConstraint('Probability_CNS').drop;

        ampl.getConstraint('ObsEq_Gap_CNS_Abs_Pos').drop;
        ampl.getConstraint('ObsEq_Gap_CNS_Abs_Neg').drop;


        ampl.getConstraint('MTR_CNS').drop;
        ampl.getConstraint('ST_CNS').drop;
        ampl.getConstraint('STL_CNS').drop;
        ampl.getConstraint('STR_CNS').drop;
        ampl.getConstraint('TIV_CNS').drop;

        ampl.getConstraint('MTS_CNS').drop;
        ampl.getConstraint('MTS_Gap_Eq_CNS_Abs_Pos').drop;
        ampl.getConstraint('MTS_Gap_Eq_CNS_Abs_Neg').drop;
        ampl.getConstraint('MTS_Gap_Eq_Sample_Abs').drop;
    end
end
