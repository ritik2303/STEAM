function [x, vds, vgs] = loadBSIMSlice(discretization_step, sampling_args)
%function [x, vds, vgs] = loadBSIMSlice(discretization_step, sampling_args)
% Author: Archit Gupta (Oct 26, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script takes in some simulation parameter and returns 1D simulation data
% from the BSIM model which is required for the compressive sensing experiment.
% The model used for simulation is always BSIM3 (v3.2.4) generated from VAPP. In
% my setup, I have a function handle for this model named bsim3.
% Please rename the function handle if you have a different setup.
%
% INPUTS:
%   discretization_step: currently, a single value of discretization step is
%       used for both vgs and vgs, so a single scalar value has to be provided.
%   sampling args (optional): A struct consisting of the following subfields
%       vds_offset: mean value of vds. The disretization for vds is done in -+
%           3V of this offset value
%       vgs_offset: mean value of vgs. The discretization for bgs is done in -+
%           3V of this offset value
%       VMIN (optional): lower boundary of discretization region for vds or vgs.
%           default -3V
%       VMAX (optional): lower boundary of discretization region for vds or vgs.
%           default 3V
%       l (optional): To be supplied when a 1-dimensional slice in the vds-vgs
%           jplane has to be sampled (scaling factor for vds)   
%       m (optional): To be supplied when a 1-dimensional slice in the vds-vgs
%           plane has to be sampled (scaling factor for vgs)   
%
%   The values of l and m determine the slice in the vds-vgs plane along which
%   the model must be evaluated. This slice is given by
%
%                       l*vds + m*vgs = 1
%
% OUTPUTS:
%   x: The data vector to be experimented with
%   vds: data points in the vds dimension
%   vgs: data points in the vgs dimension
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global STEAM_DATA_DIR;          % This variable is set by start_table_MAPP 
    if (~STEAM_DATA_DIR)
        fprintf(2, ['WARNING: Failed to read the data directory for STEAM\n.', ...
            'Setting default value!\n']);
    end

    switch(nargin)
    case 0
        discretization_step = 0.01;
        v2struct(BSIMSamplingArgs());
    case 1
        v2struct(BSIMSamplingArgs());
    case 2
        v2struct(sampling_args);
    end

    % MOD = BSIM3v3_2_4_ModSpec();    % MODSPEC model originaly in MAPP
    MOD = setMOSFETParms(bsim3('CS'), 'VA_BSIM3_NMOS');

    k = [VMIN:discretization_step:VMAX]';    % VDS range
    nSamples = length(k);

    % Sampling along the line m*vgs + l*vds = 1
    if (m ~= 0)
        vds = k + vds_offset;
        vgs = vgs_offset + (1 - l*vds)/m;
    elseif (l ~= 0)
        vgs = k + vgs_offset;
        vds = vds_offset + (1 - m*vgs)/l;
    else
        fprintf(2, ['One-dimesional slice required along lx + my = 0.\n', ...
            'However, both l, m have been supplied 0. Aborting ...\n']);
        return;
    end
    base_filename = strcat('CS_BSIM_1D_expt_lm_', num2str(l), num2str(m), ...
        '_u', num2str(1000*discretization_step), 'mV.mat');     
    data_filename = strcat(STEAM_DATA_DIR, base_filename);
    
    if (~exist(data_filename, 'file'));
        nExplicitOuts = length(MOD.ExplicitOutputNames(MOD));
        nImplicitEqns = length(MOD.ImplicitEquationNames(MOD));
        fi = zeros(nSamples, nImplicitEqns);
        qi = zeros(nSamples, nImplicitEqns);

        flag.fe = 0; flag.qe = 0; flag.fi = 1; flag.qi = 1;
        fprintf(2, 'Evaluating Model on linear grid on size %d.\n', nSamples);
        for k_index = 1:nSamples
            [~, ~, fi_k, qi_k] = MOD.fqei([vds(k_index); vgs(k_index); 0], ...
                [vds(k_index); 0], [], flag, MOD);
            fi(k_index,:) = fi_k';
            qi(k_index,:) = qi_k';
        end
        save(data_filename, 'fi', 'qi');
    else
        fprintf(2, 'Found model data in %s.\nSkipping data generation.\n', ...
            base_filename);
        load(data_filename);
    end

    x = fi(:,1); % Ids for BSIM
end