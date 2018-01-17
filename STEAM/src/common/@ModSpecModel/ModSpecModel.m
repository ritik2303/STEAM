classdef ModSpecModel < handle
    properties (SetAccess = protected, GetAccess = public)
        identifier = '';
        parm_names = {};
        parm_vals  = {};
    end
    
    properties (SetAccess = protected, GetAccess = public)
        % Variable Names
        implicit_eqn_names    = cell(0);
        explicit_output_names = cell(0);
        internal_unk_names    = cell(0);
        other_IO_names        = cell(0);
        IO_names              = cell(0);

        % Number of variables (A little redundant but this is used more often
        % than the actual names of the unknowns)
        n_explicit_outs       = 0;
        n_implicit_eqns       = 0;
        n_other_ios           = 0;
        n_internal_unks       = 0;
        support_initlimiting  = 0;

        % Other data structures
        NIL   = struct();
        J_def = struct();
        flag  = struct();
    end

    methods (Access = protected)    % Internal functions
        [fqout, dfqout] = fqei_dfqeiXYU(MOD, x_in, forq,  eori,  XorY);
        out             = functionNotDefined(MOD, eori);
    end

    methods (Abstract, Access = public)
    % ModSpecModel class defines an API that all the models will have to satisfy
    % to work with STEAM.
    %   - Currently, the functions that are supported are fqei and fqeiJ. 
    %   - Other API functions (fe, qe, ..., dfe_dvecX etc. are all performed by
    %     calling either fqei or fqeiJ in the background)
        [fqei_vals, J]   = fqeiJ(mod, x_in);
        [fe, qe, fi, qi] = fqei(mod, x_in);
    end

    methods (Access = public)
        fe = fe(MOD, VecX, VecY, VecU, ~);
        qe = qe(MOD, VecX, VecY, VecU, ~);
        fi = fi(MOD, VecX, VecY, VecU, ~);
        qi = qi(MOD, VecX, VecY, VecU, ~);
        dfe_dvecX = dfe_dvecX(MOD, VecX, VecY, VecU, ~);
        dfe_dvecY = dfe_dvecY(MOD, VecX, VecY, VecU, ~);
        dfe_dvecU = dfe_dvecU(MOD, VecX, VecY, VecU, ~);
        dqe_dvecX = dqe_dvecX(MOD, VecX, VecY, ~);
        dqe_dvecY = dqe_dvecY(MOD, VecX, VecY, ~);
        dfi_dvecX = dfi_dvecX(MOD, VecX, VecY, VecU, ~);
        dfi_dvecY = dfi_dvecY(MOD, VecX, VecY, VecU, ~);
        dfi_dvecU = dfi_dvecU(MOD, VecX, VecY, VecU, ~);
        dqi_dvecX = dqi_dvecX(MOD, VecX, VecY, ~);
        dqi_dvecY = dqi_dvecY(MOD, VecX, VecY, ~);
        % Derivative of fqei w.r.t the limited variable? This seems to make a little
        % more sense than the original derivative. Can we apply the chain rule?
        %                       vecLim = g(vecX, vecY); 
        dfe_dvecLim = dfe_dvecLim(MOD, VecX, VecY, VecU, ~);
        dqe_dvecLim = dqe_dvecLim(MOD, VecX, VecY, VecU, ~);
        dfi_dvecLim = dfi_dvecLim(MOD, VecX, VecY, VecU, ~);
        dqi_dvecLim = dqi_dvecLim(MOD, VecX, VecY, VecU, ~);

        % TODO: Please remove this function if it is not needed. Run tests after
        % commenting the next line to confirm.
        dfqei_dvecX = dfqei_dvecX(MOD, VecX, VecY, VecU, ~);

        % Getting and Setting the parameters
        parms = getparms(MOD, parm_name, ~);
        oMOD  = setparms(MOD, parm_names, parm_vals, ~);

        % Cannot ignore limiting for any of the cases. It is pretty much
        % essential to get any device working Look up DAAV6ModSpec.m for
        % reference on how to write these functions on your own. The
        % ModSpec_common_skeleton defines these functions using the vecvalder
        % object which is the source of all the problems. So, we are going to
        % redefine them all
    
        dlimiting_dvecX = dlimiting_dvecX(MOD, VecX, VecY, VecU, ~);
        dlimiting_dvecY = dlimiting_dvecY(MOD, VecX, VecY, VecU, ~);

        function mod = ModSpecModel(o_model)
        % Instantiating the ModSpecModel class from an existing model in MAPP
            if ( isstr(o_model) )
                mod.identifier = o_model;
            elseif ( isa(o_model, 'ModSpecModel') )
                % [TODO]: Write an explicit copy method and use it
                mod.identifier            = o_model.identifier;

                % Variable names
                mod.explicit_output_names = o_model.explicit_output_names;
                mod.other_IO_names        = o_model.other_IO_names;
                mod.IO_names              = o_model.IO_names;

                % Number of variables of each different type
                mod.n_explicit_outs = o_model.n_explicit_outs;
                mod.n_implicit_eqns = o_model.n_implicit_eqns;
                mod.n_other_ios     = o_model.n_other_ios;

                % Init/Limiting support
                mod.support_initlimiting = o_model.support_initlimiting;

                % Parameters
                mod.parm_names = o_model.parm_names;
                mod.parm_vals  = o_model.parm_vals;

                % Network Interface Layer
                mod.NIL = o_model.NIL;
            else
                % Assume that a MODSPEC Model has been supplied
                mod.identifier = o_model.uniqID;
                mod.parm_names = o_model.parm_names;
                mod.parm_vals  = o_model.parm_vals;

                % Number of variables of each different type
                mod.n_implicit_eqns = length( o_model.ImplicitEquationNames(o_model) );
                mod.n_explicit_outs = length( o_model.ExplicitOutputNames(o_model) );
                mod.n_other_ios     = length( o_model.OtherIONames(o_model) );
                mod.n_internal_unks = length( o_model.InternalUnkNames(o_model) );

                % Variable names
                mod.other_IO_names = o_model.OtherIONames(o_model);
                mod.explicit_output_names = o_model.ExplicitOutputNames(o_model);
                mod.internal_unk_names = o_model.InternalUnkNames(o_model);

                % [TODO] This is certainly not enough. Do something with the
                % model. Atleast copy over the important details / function
                % handles
            end 
            mod.J_def = FQJAC(mod);
            mod.flag  = Flag();
        end

        function oios = OtherIONames(MOD, ~)
            oios = MOD.other_IO_names;
        end 

        function eos = ExplicitOutputNames(MOD, ~)
            eos = MOD.explicit_output_names;
        end

        function ios = IOnames(MOD, ~)
            ios = { MOD.other_IO_names{:}, MOD.explicit_output_names{:} };
        end

        function iunks = InternalUnkNames(MOD, ~)
            iunks = MOD.internal_unk_names;
        end

        function ieqs = ImplicitEquationNames(MOD, ~)
            ieqs = MOD.implicit_eqn_names;
        end

        function unames = uNames(MOD, ~)
            unames = {};
        end

        function pnames = parmnames(MOD, ~) % What's with the naming convention... We had camel-case so far
            pnames = MOD.parm_names;
        end

        function pdefs = parmdefaults(MOD, ~)
            % This is a tough one: HACK [TODO]: Just return parm_vals
            pdefs = MOD.parm_vals;
        end
    end
end
