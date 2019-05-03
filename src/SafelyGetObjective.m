%*******************************************************************************
% SafelyGetObjective
%
% 05/01/19:
% Not sure if this is still relevant with new versions of the AMPL-Matlab API,
% but it is harmless, so might as well leave it in
%
% (Some time long ago):
% Very annoyingly, I've found that on some rare occasions one can get a strange
% Java error when trying to pull the value of an objective function from AMPL.
% Even stranger, a solution seems to be to create an ampl instance of the
% objective and then try to display it before pulling the value.
% No idea what's going on and haven't been able to construct a MWE.
% But this bit of code seems to take care of the problem.
%*******************************************************************************
function val = SafelyGetObjective(ampl, ObjectiveName)
    try
        val = ampl.getValue(ObjectiveName);
    catch
        Obj = ampl.getObjective(ObjectiveName);
        evalc('ampl.display(Obj)');

        val = ampl.getValue(ObjectiveName);
    end
end
