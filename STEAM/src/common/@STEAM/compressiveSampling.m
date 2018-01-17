function lookUpTable = compressiveSampling(sweepSizes, samplePoints, A, ...
    flag, MOD)
% function lookUpTable = compressiveSampling(sweepSizes, samplePoints, A, ...
%     flag, MOD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 SAMPLING AND EVALUATING THE DEVICE MODEL                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a randomly sampled set of points, this function samples the MOSFET at
% these points and calls CS-based reconstruction algorithm to fill in the
% the values at the entire rectangular gird
    nSamplePoints = length(samplePoints{1});
    nTabulatedVars = MOD.nTabulatedVars;
    nInternalUnks = MOD.nInternalUnks;
    nExplicitVars = min(2,MOD.nExplicitVars);
    nOtherIOVars  = min(2, MOD.nOtherIOVars);
    nImplicitEqns = length(MOD.ImplicitEquationNames(MOD));

    sampledTable = struct();
    if (flag.fe)
        sampledTable.fe = zeros(nSamplePoints, nExplicitVars);
        sampledTable.qe = zeros(nSamplePoints, nExplicitVars);
        if (flag.fi)
            % This will not be used for  MOSFETs. This part of the code was
            % written to handle 2-terminal devices with internal unknowns. Since
            % we do not support 3 or higher dimensional splines at the moment,
            % this is required to get some variety of devices to work.
            sampledTable.fi = zeros(nSamplePoints, nImplicitEqns);
            sampledTable.qi = zeros(nSamplePoints, nImplicitEqns);
        end
    else
        fprintf(2,'WARNING: No explicit variables for tabulation\n');
    end

    run4TMOSFETHacks = flag.fi;

    if (nOtherIOVars == 2);
        if (run4TMOSFETHacks)
            for sample_index = 1:nSamplePoints
                VecX = [samplePoints{1}(sample_index); samplePoints{2}(...
                    sample_index); 0];
                VecY = [VecX(1,1); 0];
                [fe, qe, fi, qi] = MOD.base_fqei(VecX, VecY, [], [], ...
                    flag, MOD);
                sampledTable.fe(sample_index,:) = [fe(1,1)+fi(1,1); fe(2,1)];
                sampledTable.qe(sample_index,:) = [qi(1,1);qe(2,1)]; 
            end 
            flag.fi = 0;    % We don't need the flag anymore, the new model (and
            % the following compressed sensing algorithm will only work on fe
        else
            fprintf(2, ['You have managed to drag yourself into some very', ...
                'deep mess. Kindly contact the author for details\n']);
        end
    else
        fprintf(2, ['Sorry, compressive sampling is only available for ', ...
            'MOSFETs at the moment\n']);
        return;
    end

    grid_shaped_feqe_dims = cat(2, sweepSizes, [nExplicitVars]);
    fe_CS = zeros(grid_shaped_feqe_dims);
    qe_CS = zeros(grid_shaped_feqe_dims);
    if (flag.fe)    % Reconstruct fe, qe
        data_dims = size(samplePoints{1});
        for recon_dim = 1:nExplicitVars
            fe_CS(:,:,recon_dim) = reshape(reconstruct_with_CS(sampledTable.fe(:, ...
                recon_dim), A, sweepSizes), sweepSizes);
            qe_CS(:,:,recon_dim) = reshape(reconstruct_with_CS(sampledTable.qe(:, ...
                recon_dim), A, sweepSizes), sweepSizes);
        end
        lookUpTable.fe = fe_CS;
        lookUpTable.feqe = cat(3, fe_CS, qe_CS);
        if (flag.fi)
            fprintf(2, ['Reconstruction of fi for models is not supported', ...
                'with CS at the moment\n']);
        end
    end
end
