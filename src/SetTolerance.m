%*******************************************************************************
% SetTolerance
%
% Adjust solver tolerance without disturbing the rest of the settings.
% This solver tolerance is only for CNS
%*******************************************************************************
function SetTolerance(ampl, Tolerance)

    Solver = ampl.getOption('solver');

    if strcmp(Solver, 'cplex')
        SolverOptionsName = 'cplex_options';
        OptionName = 'feasibility';
    elseif strcmp(Solver, 'gurobi')
        SolverOptionsName = 'gurobi_options';
        OptionName = 'feastol';
    else
        ErrMsg = sprintf('No code for SetTolerance for solver %s', Solver);
        error(ErrMsg);
    end

    % ampl.getOption is a Java string, so convert it
    Options = char(ampl.getOption(SolverOptionsName));
    % [^\s]+ means find everything up to a whitespace character
    Find = [OptionName '=[^\s]+'];
    Replace = [OptionName '=' num2str(Tolerance)];
    if isempty(regexp(Options, Find))
        NewOptions = [Options ' ' Replace];
    else
        NewOptions = regexprep(Options, Find, Replace);
    end
    ampl.setOption(SolverOptionsName, NewOptions)
end
