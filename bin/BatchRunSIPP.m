function [] = BatchRunSIPP(SaveDir, SimSet)

    errstr = 'Need to pass nonempty SaveDir for this routine.';
    if ~exist('SaveDir', 'var')
        error(errstr);
    else
        if isempty(SaveDir)
            error(errstr);
        end
    end

    sbase = ['!matlab -nodesktop -nosplash -singleCompThread'...
         ' -r "RunSIPP(''%s'', ''%s'', %d, 1)" &'];

    if ~exist('SimSet')
        error('Must pass SimSet')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SOME HARD-CODING AHEAD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NSIMSROBUST = 9;
    switch SimSet
        case 'main'
            NSIMS = 12;
        case 'sigma'
            NSIMS = NSIMSROBUST;
        case 'sigma-young'
            NSIMS = NSIMSROBUST;
        case 'extra'
            NSIMS = 8;
        otherwise
            error('Case not recognized.')
    end

    for i = 1:1:NSIMS
        s = sprintf(sbase, SaveDir, SimSet, i);
        disp(s);
        eval(s);
        pause(10);
    end
end
