function subcktdata = MOSInverter(model, parm_string, subCkt, method)
% function subcktdata = MOSInverter(model, parm_string, subCkt, method)

    subcktdata.cktname = 'MOS Inverter';
    do_tabulate = 0;

    if (nargin > 3)
        do_tabulate = 1;
        subcktdata.cktname = 'Tabulated MOS Inverter';
    end

    pmos_parm_string = strcat(parm_string, '_PMOS');
    nmos_parm_string = strcat(parm_string, '_NMOS');

    % nodes (names)
    subcktdata.nodenames = {'vdd', 'in', 'out'};
    subcktdata.terminalnames = subcktdata.nodenames;  % Used in ring-oscillator
    subcktdata.groundnodename = 'gnd';
    subcktdata.elements = {};

    % list of elements 

    % pmosElem
    if (do_tabulate)
        subcktdata = add_subcircuit(subcktdata, subCkt(model, ...
            pmos_parm_string, method), 'PMOS', {'out', 'in', 'vdd'});
    else
        subcktdata = add_subcircuit(subcktdata, subCkt(model, ...
            pmos_parm_string), 'PMOS', {'out', 'in', 'vdd'});
    end

    % nmosElem
    if (do_tabulate)
        subcktdata = add_subcircuit(subcktdata, subCkt(model, ...
            nmos_parm_string, method), 'NMOS', {'out', 'in', 'gnd'});
    else
        subcktdata = add_subcircuit(subcktdata, subCkt(model, ...
            nmos_parm_string), 'NMOS', {'out', 'in', 'gnd'});
    end
end
