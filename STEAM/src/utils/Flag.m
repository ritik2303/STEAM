classdef Flag < handle
    properties (SetAccess = private, GetAccess = public)
        fe = 0;
        fi = 0;
        qe = 0;
        qi = 1;
    end

    properties (Access = public) % Look into fqei_from_fqeiJ_ModSpec for reasons - AG
        J = 0;
    end

    methods (Access = public)
        function Obj = Flag()
            Obj.setFQ();
            Obj.J = 0;
        end

        function setJ(Obj)
            Obj.J = 1;
        end

        function clearJ(Obj)
            Obj.J = 0;
        end

        function clearFQ(Obj)
            Obj.fe = 0;
            Obj.qe = 0;
            Obj.qi = 0;
            Obj.fi = 0;
        end

        function setFQ(Obj)
            Obj.fe = 1;
            Obj.qe = 1;
            Obj.qi = 1;
            Obj.fi = 1;
        end

        function clearF(Obj)
            Obj.fe = 0;
            Obj.fi = 0;
        end

        function clearQ(Obj)
            Obj.qe = 0;
            Obj.qi = 0;
        end

        function setF(Obj)
            Obj.fe = 1;
            Obj.fi = 1;
        end

        function setQ(Obj)
            Obj.qe = 1;
            Obj.qi = 1;
        end
    end
end

