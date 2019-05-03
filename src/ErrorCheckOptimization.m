%*******************************************************************************
% ErrorCheckOptimization
%
% Ensure that the return code on the optimization was acceptable.
%*******************************************************************************
function [] = ErrorCheckOptimization(ampl, Identifier, InfeasibleFlag)
    SolveResult = ampl.getValue('solve_result');
    if (      (strcmp(SolveResult, 'solved') == 0) ...
            & (strcmp(SolveResult, 'infeasible') == 0))
        disp(sprintf(...
            ['\t Return code other than solved or infeasible.\n'...
             '\t In %s:'],...
            Identifier));
        disp(sprintf('\t solve_result = %s.', SolveResult));
        disp(sprintf('\t solve_result_num = %d.',...
                ampl.getValue('solve_result_num')));
    end
    diary off; diary on; % Flush

    if (      (strcmp(SolveResult, 'infeasible') == 1) ...
            & (InfeasibleFlag == 1))
        disp(sprintf('\t INFEASIBILITY IN %s.', Identifier));
        disp(sprintf('\t solve_result = %s.', SolveResult));
        disp(sprintf('\t solve_result_num = %d.',...
                ampl.getValue('solve_result_num')));
        disp(sprintf('\t CriterionHat = %f.', ...
                ampl.getValue('CriterionHat')));
    end
end
