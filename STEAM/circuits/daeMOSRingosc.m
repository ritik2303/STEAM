function [dae, outputs, sim_args] = daeMOSRingosc(n_elements, model, ...
    parm_string, subCkt, method)
% function [dae, outputs, sim_args] = daeMOSRingosc(n_elements, model, ...
%     parm_string, subCkt, method)
% An 'N' stage ring oscillator with BSIM/MVS MOSFET
% Author: J. Roychowdhury, 2012/07/18
% modified: Archit Gupta, 2016/03/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The Circuit: 
%	An 'N' stage ring oscillator made with MVS/BSIM MOSFETs using 0.18u
%	parameters, probably from MOSIS. The value of 'N' can be set using
%	'n_elements' input argument
% Modification:
%       Minor changes for using table-based MOSFET models instead of using the
%       regular model for simulation. Also, the original DAE had a fixed number
%       of elements (3). This has been changed to an input argument
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% see STEAM.m for a description of functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Type "help MAPPlicense" at the MATLAB/Octave prompt to see the license      %
%% for this software.                                                          %
%% Copyright (C) 2008-2013 Jaijeet Roychowdhury <jr@berkeley.edu>. All rights  %
%% reserved.                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cktdata.cktname = ['MOS Ring Oscillator (', num2str(n_elements), '-Stage)'];
    do_tabulate = 0;

    if (nargin > 4)
        do_tabulate = 1;
        cktdata.cktname = ['Tabulated ', cktdata.cktname];
    end

    vdd_dc = 1.0;
    c_l = 1e-19;

    %===========================================================================
    % subcircuit: MOSinverter
    %---------------------------------------------------------------------------

    inverter = [];
    if (do_tabulate)
        inverter = MOSInverter(model, parm_string, subCkt, method);
    else
        inverter = MOSInverter(model, parm_string, subCkt);
    end

    inverter = add_element(inverter, capModSpec(), 'c_l', ...
       {'out', 'gnd'}, {{'C', c_l}}, {});

    %===========================================================================

    cktdata.groundnodename = 'gnd';

    cktdata.nodenames = cell(1, n_elements+1);
    cktdata.nodenames{1} = 'vdd';
    for n_inv =  1:n_elements 
        % Inverter i
        cktdata.nodenames{1+n_inv} = ['inv', num2str(n_inv)];
    end

    cktdata = add_element(cktdata, vsrcModSpec(), 'Vdd', ...
    {'vdd', 'gnd'}, {}, {{'E', {'DC', vdd_dc}}});

    for n_inv = 1:n_elements-1
        % X_{i} Elem
        cktdata = add_subcircuit(cktdata, inverter, ['X', num2str(n_inv)], ...
            {'vdd', cktdata.nodenames{1+n_inv}, cktdata.nodenames{2+n_inv}});
    end
    cktdata = add_subcircuit(cktdata, inverter, ['X', num2str(n_elements)], ...
        {'vdd', cktdata.nodenames{1+n_elements}, 'inv1'});

    
    dae = MNA_EqnEngine(cktdata); 
    outputs = StateOutputs(dae);
    outputs = outputs.DeleteAll(outputs);
    outputs = outputs.Add({'e_inv1', 'e_inv2', 'e_inv3'}, outputs);

    % Finding if the model is BSIM or MVS and setting the appropriate simulation
    % arguments (The two circuits have different values of tstop/tstep etc.
    sim_args = struct();
    sim_args.tstart = 0;

    if (findstr(parm_string, 'BSIM'))
        sim_args.tstep = 6e-9;
        sim_args.tstop = 4e-7;
    elseif(findstr(parm_string, 'MVS'))
        sim_args.tstep = 6e-13;
        sim_args.tstop = 4e-11;
    elseif (findstr(parm_string, 'PSP'))
        sim_args.tstep = 1e-9;
        sim_args.tstop = 4e-7;
    else
        error('ERROR: circuit parameters not found for the given mos model\n');
    end
    sim_args.VDD = vdd_dc;
    vdd_args.VDD = vdd_dc;

    % Transient function descriptions
    vddFunc = @(t, a) a.VDD;

    % Setting up transient functions
    dae = dae.set_utransient('Vdd:::E', vddFunc, vdd_args, dae);

    % Setting up DC inputs
    QSS = dae.utransient(sim_args.tstart, dae);
    dae = dae.set_uQSS(QSS, dae);

    % Setting up inputs for HB
    dae = dae.set_uHB('Vdd:::E', vddFunc, vdd_args, dae);
end
