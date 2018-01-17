classdef FQEI < MODDIMS
    properties ( Access = public )
        fe = [];
        fi = [];
        qe = [];
        qi = []
    end

    methods ( Access = public )
        function fqei = FQEI(MOD, fe, qe, fi, qi)   % Pass in a ModSpec model to gather data
            fqei = fqei@MODDIMS(MOD);

            if ( nargin == 1)
                % Setting up Fe/Qe/Fi/Qi
                fqei.fe = zeros(fqei.n_explicit_outs, 1);
                fqei.fi = zeros(fqei.n_implicit_eqns, 1);
                fqei.qe = zeros(fqei.n_explicit_outs, 1);
                fqei.qi = zeros(fqei.n_implicit_eqns, 1);
            elseif ( nargin == 5 ) 
                % Everything has to be supplied
                fqei.fe = fe;
                fqei.fi = fi;
                fqei.qe = qe;
                fqei.qi = qi;
            else
                error( 'ERROR: Support either 1 (MODEL) or 5 (MODEL, VALUES) only\n');
            end
        end
    end
end
