function [dae, outputs, simArgs] = daeMOSInverter(varargin)
%function [dae, outputs, simArgs] = daeMOSInverter(varargin)

    cktdata.cktname = 'MOS Inverter';
    if (nargin > 2)
        cktdata.cktname = 'Tabulated MOS Inverter';
    end

    % nodes (names)
    cktdata.nodenames = {'vdd', 'in', 'out'};
    cktdata.groundnodename = 'gnd';

    % list of elements 

    % vsupElem
    vsupElem.name = 'vdd';
    vsupElem.model = vsrcModSpec('vdd');
    vsupElem.nodes = {'vdd', 'gnd'}; % p, n
    vsupElem.parms = {}; % vsrc/isrc have no parameters

    cktdata.elements = {vsupElem};

    % vinElem
    vinElem.name = 'vin';
    vinElem.model = vsrcModSpec('vin');
    vinElem.nodes = {'in', 'gnd'}; % p, n
    vinElem.parms = {}; % vsrc/isrc have no parameters

    cktdata.elements = {cktdata.elements{:}, vinElem};

    cktdata = add_subcircuit(cktdata, MOSInverter(varargin{:}), 'INV', ...
        {'vdd', 'in' , 'out'});
    
    dae = MNA_EqnEngine(cktdata);

    % Setting up transient/DC/AC parameters

    VDD = 1.0;
    VInOffset = 0.5;

    % Finding if the model is BSIM or MVS
    simArgs = struct();
    simArgs.fstart = 1;
    simArgs.fstop = 1e9;
    simArgs.nsteps = 5; % 5 points per decade

    parm_string = varargin{2};
    kBSIM = findstr(parm_string, 'BSIM');
    kMVS = findstr(parm_string, 'MVS');
    
    % AC analysis arguments
    simArgs.VIN_MIN = VDD/2 - 0.1;
    simArgs.VIN_MAX = VDD/2 + 0.1;
    simArgs.n_ops = 4;
    
    % Different Transient Inputs to try (Declare here as they are required for setting the time step)
    vddArgs.VDD = VDD;

    sinArgs.A = .2;
    sinArgs.f = 1e6;
    sinArgs.Offset = VInOffset;

    if (~isempty(kBSIM))
        simArgs.tstop = 4/sinArgs.f;

        % DC Analysis Arguments
        simArgs.v_start = 0;
        simArgs.v_step = 0.01;
        simArgs.v_stop = VDD;

        % Transient Arguments
        simArgs.tstart = 0;
        simArgs.tstep = 0.05/sinArgs.f;
        simArgs.tstop = 4/sinArgs.f;
    elseif(~isempty(kMVS))
        sinArgs.A = 0.2;

        % Transient Arguments
        simArgs.tstart = 0;
        simArgs.tstep = 0.05/sinArgs.f;
        simArgs.tstop = 4/sinArgs.f;

        % DC Analysis Arguments
        %simArgs.v_start = 0.85*VDD;
        simArgs.v_start = VDD;
        simArgs.v_step = -0.01;
        simArgs.v_stop = 0;
    else
        fprintf(2, 'Circuit parameters not found for the given MOS model\n');
    end

    % This transient input requires variables that were described above
    stepArgs.A = .2;
    stepArgs.scale = 10*simArgs.tstep;
    stepArgs.Offset = 0.4;
    stepArgs.Onset = (simArgs.tstart+simArgs.tstop)/2;
    stepArgs.smoothing = 0.01;

    % Transient function description
    sinFunc = @(t,args) args.Offset+args.A*sin(args.f*2*pi*t);
    stepFunc = @(t,args) args.Offset+args.A*smoothstep((t-args.Onset)/...
        args.scale, args.smoothing);

    vddFunc = @(t,a) a.VDD;
    
    inFunc = sinFunc; inArgs = sinArgs;
    %inFunc = stepFunc; inArgs = stepArgs;

    % Setting up transient inputs
    dae = dae.set_utransient('vin:::E', inFunc, inArgs, dae);
    dae = dae.set_utransient('vdd:::E', vddFunc, vddArgs, dae);
    
    % Setting up DC inputs
    QSS = dae.utransient(simArgs.tstart, dae);
    dae = dae.set_uQSS(QSS, dae); % Vin is irrelevant for DCSWEEP

    outputs = StateOutputs(dae);
    outputs = outputs.DeleteAll(outputs);
    outputs = outputs.Add({'e_in', 'e_out'}, outputs);
    % outputs = outputs.Add({'e_out'}, outputs); % For AC Analysis
end
