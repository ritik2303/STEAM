function subcktnetlist = MOSFET3T_subckt(model, parm_string, method)
% This file is truly needed only to make the example scripts easier to use.
% It allows all MOSFET models to be usable with the same calling syntax. To
% do this, we make a subcircuit even out of individual 3-Terminal MOS models.
    subcktnetlist.cktname = 'MOSFET_3T';
    subcktnetlist.nodenames = {'d', 'g', 's'};
    subcktnetlist.groundnodename = 'gnd';

    % Adding the core MOS Model
    MOSElem.name = parm_string;
    if (nargin > 2)
        MOSElem.model = STEAM(MOSModel2D(model, parm_string), method.i_method, ...
            method.d_method, method.order, method.bounds);
    else
        MOSElem.model = MOSModel2D(model, parm_string);
    end
    MOSElem.nodes = {'d', 'g', 's'};
    MOSElem.parms = MOSElem.model.getparms();
    subcktnetlist.elements = {MOSElem};

    subcktnetlist.terminalnames = {'d', 'g', 's'};
end