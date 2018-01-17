classdef NIL < handle
    properties (SetAccess = private, GetAccess = public)
        node_names   = cell(0);
        refnode_name = '';
        io_types     = cell(0);
        io_nodenames = cell(0);
    end

    methods (Access = public)
        function Obj = NIL(node_names)
        % Class constructor
            Obj.node_names   = node_names;
            Obj.refnode_name = node_names{end};
            Obj.io_types     = cell(1, 2 * (length(Obj.node_names) - 1));
            Obj.io_nodenames = cell(1, 2 * (length(Obj.node_names) - 1));

            % Currently, with STEAM, we support only electrical devices
            % Set up Input/Output types to 'voltages', and 'currents'
            [Obj.io_types{1:length(Obj.node_names)-1}] = deal('v');
            [Obj.io_types{length(Obj.node_names):end}] = deal('i');

            % Names of nodes (auto-generated)
            [Obj.io_nodenames{1:2:end}] = deal(Obj.node_names{1});
            [Obj.io_nodenames{2:2:end}] = deal(Obj.node_names{2});
        end

        function n_names = NodeNames(Obj, ~)
            n_names = Obj.node_names;
        end

        function r_name = RefNodeName(Obj, ~)
            r_name = Obj.refnode_name;
        end

        function io_types = IOtypes(Obj, ~)
            io_types = Obj.io_types;
        end

        function io_nodenames = IOnodeNames(Obj, ~) % Why is the first n not capital here - AG
            io_nodenames = Obj.io_nodenames;
        end
    end
end
