function [x, vds, vgs] = loadBSIMImage(discretization_step, sampling_args)
%function [x, vds, vgs] = loadBSIMImage(discretization_step, sampling_args)
% Author: Archit Gupta (Oct 26, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script takes in some simulation parameter and returns 2D simulation data
% from the BSIM model which is required for the compressive sensing experiment.
% The model used for simulation is always BSIM3 (v3.2.4) generated from VAPP. In
% my setup, I have a function handle for this model named bsim3.
% Please rename the function handle if you have a different setup.
%
% INPUTS:
%   discretization_step: currently, a single value of discretization step is
%   used for both vgs and vgs, so a single scalar value has to be provided.
%   sampling args (optional): A struct consisting of the following subfields
%       vds_offset: mean value of vds. The disretization for vds is done in -+ 3V
%           of this offset value
%       vgs_offset: mean value of vgs. The discretization for bgs is done in -+
%           3V of this offset value
%       VMIN (optional): lower boundary of discretization region for vds or vgs.
%           default -3V
%       VMAX (optional): lower boundary of discretization region for vds or vgs.
%           default 3V
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global STEAM_DATA_DIR;          % This variable is set by start_table_MAPP 
    if (~STEAM_DATA_DIR)
        fprintf(2, ['WARNING: Failed to read the data directory for STEAM\n.', ...
            'Setting default value!\n']);
    end

    default_sampling_args = BSIMSamplingArgs();

    switch(nargin)
    case 0
        discretization_step = 0.1;
        v2struct(default_sampling_args);
    case 1
        v2struct(default_sampling_args);
    case 2
        v2struct(sampling_args);
    end

    % MOD = BSIM3v3_2_4_ModSpec();    % MODSPEC model originaly in MAPP
    MOD = setMOSFETParms(bsim3('CS'), 'VA_BSIM3_NMOS');

    k = [VMIN:discretization_step:VMAX]';    % VDS range
    n_samples = length(k);

    vds = vds_offset + k;
    vgs = vgs_offset + k;
    image_dimensions = {length(vds) length(vgs)};

    base_filename = strcat('CS_BSIM_2D_expt_lm_', num2str(l), num2str(m), ...
        '_u', num2str(1000*discretization_step), 'mV.mat');     
    data_filename = strcat(STEAM_DATA_DIR, base_filename);
    
    if (~exist(data_filename, 'file'));
        n_explicit_outs = length(MOD.ExplicitOutputNames(MOD));
        n_implicit_eqns = length(MOD.ImplicitEquationNames(MOD));
        fi = zeros(image_dimensions{:}, n_implicit_eqns);
        qi = zeros(image_dimensions{:}, n_implicit_eqns);

        flag.fe = 0; flag.qe = 0; flag.fi = 1; flag.qi = 1;
        fprintf(2, 'Evaluating Model on linear grid on size %dx%d\n', n_samples, n_samples);
        for vds_index = 1:length(vds)
            for vgs_index = 1:length(vgs)
                [~, ~, fi_k, qi_k] = MOD.fqei([vds(vds_index); ...
                    vgs(vgs_index); 0], [vds(vds_index); 0], [], flag, MOD);
                fi(vds_index, vgs_index, :) = fi_k';
                qi(vds_index, vgs_index, :) = qi_k';
            end
        end

        save(data_filename, 'fi', 'qi');
    else
        fprintf(2, 'Found model data in %s.\nSkipping data generation\n', ...
            base_filename);
        load(data_filename);
    end

    x = squeeze(fi(:,:,1));
end