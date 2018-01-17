function lookUpTable = uniformSampling(sweepSizes, inputSweep, flag, MOD)

    nTabulatedVars = MOD.nTabulatedVars;
    nInternalUnks = MOD.nInternalUnks;
    nExplicitVars = min(2,MOD.nExplicitVars);
    nOtherIOVars  = min(2,MOD.nOtherIOVars);
    nImplicitEqns = length(MOD.ImplicitEquationNames(MOD));

    run4TMOSFETHacks = 1;   % We can definitely obtain the values of fe/fi by
                            % running a DC sweep
    % Initializing and setting up the look-up table itself
    lookUpTable = struct();
    if (flag.fi && ~run4TMOSFETHacks)
        lookUpTable.fefi = zeros([sweepSizes nExplicitVars+nImplicitEqns]);
        lookUpTable.fqei = zeros([sweepSizes 2*(nExplicitVars+nImplicitEqns)]);
    elseif (flag.fe)
        lookUpTable.fe = zeros([sweepSizes nExplicitVars]);
        lookUpTable.feqe = zeros([sweepSizes 2*nExplicitVars]);
    else
        fprintf(2,'WARNING: No explicit variables for tabulation\n');
    end

    % Generalizing the method for models with arbitrary number of internal
    % unknowns (but less than or equal to 2 Explicit outputs (i.e. size of fe/qe
    % functions should be <= 2). We will do this by running DC sweep over the
    % given range of VecX -- OtherIO values


    if (nOtherIOVars == 1)
        if (run4TMOSFETHacks)
            if (nInternalUnks > 1)
                fprintf(2, ['ERROR: Tabulation of 1D models with', ...
                    '#InternalUnknowns>1 not supported\n']);
                return;
            else
                for xi=1:sweepSizes(1)
                    VecX = inputSweep{1}(xi);
                    for xj=1:sweepSizes(2)
                        VecY = inputSweep{2}(xj);
                        [fe, qe, fi, qi] = MOD.base_fqei(VecX, VecY, [], [], ...
                            flag, MOD);
                        lookUpTable.fefi(xi,xj,:) = [fe; fi];
                        lookUpTable.fqei(xi,xj,:) = [fe; fi; qe; qi];
                    end
                end
            end
        else
            for xi=1:sweepSizes(1)
                [fe, qe, fi, qi] = MOD.base_fqei(inputSweep{1}(xi), [], [], ...
                     flag, MOD);
                % TODO: In this case, the only variable cannot be an Internal
                % Unknown so there is no problem, but we cannot generalize this to
                % the system where we have a total of 2 variables
                lookUpTable.fe(xi,:) = fe;
                lookUpTable.feqe(xi,:) = [fe; qe];
            end
        end
    elseif (nOtherIOVars == 2)
        if (run4TMOSFETHacks) % HACK: Specific for BSIM/MVS
            % TODO: This is a really hack version only for BSIM/MVS models. need
            % to write a better (and somewhat theoritically correct version)
            % version involving DC sweep for estimating the new fe(vecX) and
            % probably ??? Transient/AC analysis for qe(vecX). Right now, we
            % will make use of the knowledge of the model structure (Internal
            % nodes represent the nodes at which resistances are connected for
            % the resistances at drain and source 
            VecX = [0; 0; 0];   % VecX(3,1) = vsb is always 0 (source/bulk tied)
            for xi=1:sweepSizes(1)
                VecX(1,1) = inputSweep{1}(xi);
                for xj=1:sweepSizes(2)
                    VecX(2,1) = inputSweep{2}(xj);
                    [fe, qe, fi, qi] = MOD.base_fqei(VecX, [VecX(1,1); 0], ...
                        [], flag, MOD);
                    lookUpTable.fe(xi,xj,:) = [fe(1,1) + fi(1,1); fe(2,1)];  % I{d,g}s
                    % Sign change because the current going into the outer node is
                    % same as the current going out of the internal node
                    lookUpTable.feqe(xi,xj,:) = [fe(1,1)+fi(1,1); fe(2,1); ...
                        qi(1,1);qe(2,1)]; 
                end
            end
        else 
            VecX = [];
            for xi=1:sweepSizes(1)
                VecX(1,1) = inputSweep{1}(xi);
                for xj=1:sweepSizes(2)
                    VecX(2,1) = inputSweep{2}(xj);
                    [fe, qe, ~, ~] = MOD.base_fqei(VecX, [], [], flag, MOD);
                        % This part of the code deals with a model that has no
                        % internal unknowns. Here, things are simpler (and also
                        % have some sort of a firm basis). Just evaluate fe/qe
                        % and copy them over to the table. fi/qi are not in the
                        % picture
                    lookUpTable.fe(xi,xj,:) = fe;
                    lookUpTable.feqe(xi,xj,:) = [fe; qe];
                end
            end
        end
    else
      fprintf(2,'ERROR: tabulate_ModSpec does not support input size > 3\n');
      lookUpTable = [];
    end
end

