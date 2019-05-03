%*******************************************************************************
% InitializeAMPL
%
% Start an AMPL instance and load model files
%*******************************************************************************
function [ampl] = InitializeAMPL(ModelFiles, Settings)
    ampl = AMPL;

    for m = 1:1:length(ModelFiles)
        if ~exist(ModelFiles{m}, 'file')
            error('Model file %d could not be found.', m);
        else
            % The exist looks through the path--but need to translate for AMPL
            fn = which(ModelFiles{m});
            ampl.read(fn);
        end
    end

    SetAMPLOptions(ampl, Settings);
end

%*******************************************************************************
% SetAMPLOptions
%*******************************************************************************
function SetAMPLOptions(ampl, Settings)
    % Some innocuous settings that I do not allow the user to change.
    ampl.setOption('TMRDIR', '"./"');
    ampl.setOption('presolve_warnings', '-1');
    ampl.setOption('solver_msg', '1');
    ampl.setOption('print_precision', '6');
    ampl.setOption('print_round', '5');
    ampl.setOption('show_stats', '1');

    % Solver settings
    ampl.setOption('gurobi_options',...
        ['outlev=1 '...
         'timing=1 '...
                            ]);
    ampl.setOption('cplex_options', ...
        ['display=1 '...
         'parallelmode=1 '...
                                ]);

    % AMPL presolve tolerance
    if isfield(Settings, 'PreSolveEps')
        if ~isempty(Settings.PreSolveEps)
            ampl.setOption('presolve_eps',...
                num2str(Settings.PreSolveEps));
        end
    end

    % Set solver
    if isfield(Settings, 'Solver')
        if ~ismember(Settings.Solver, {'cplex', 'gurobi'})
            error('Solver needs to be either cplex or gurobi.')
        else
            ampl.setOption('solver', Settings.Solver);
        end
    end
end
