function ODEXYobj = ODEXY(DAE, TRmethod, tranparms) % DAE=DAEAPIv6.2
%function ODEXYobj = ODEXY(DAE, TRmethod, tranparms) % DAE=DAEAPIv6.2
%Author: Jaijeet Roychowdhury <jr@berkeley.edu> 2008-2013
%
%ODEXY: transient simulation object using Matlab's built-in ODE/DAE solvers.
%
%arguments:
% 	- DAE: mandatory
% 	- TRmethod: optional. Can be one of:
% 	  - 'ode15s' (or []): this is the default
%	    - warning: no NR hooks (like init/limiting) available
% 	  - 'ode45' - ODEs only
% 	  - 'ode23' - ODEs only
% 	- tranparms: structure used to set DAE/ODE solution parameters.
% 	       The following fields are supported (help odeset for
% 	       descriptions of most of these):
% 	       - tranparms.RelTol: defaults to 1e-4
% 	       - tranparms.AbsTol: defaults to 1e-9
% 	       - tranparms.MaxStep: default is (tstop-tstart)/10
% 	       - tranparms.InitialStep: defaults to [].
% 	       - tranparms.InitialSlope: defaults to a zero vector.
% 	       - tranparms.BDF: ['on' | 'off' ] (default is 'off')
% 	       - tranparms.MaxOrder: [ 1|2|3|4|5 ] (default is 5)
% 	       - tranparms.BEforFirstStep: [ 0|1 ] (default is 1)
%		 - use Backward Euler for the first timestep (hack to
%		   make DAE initial conditions consistent - doesn't
%		   always help).
%	         - note: this overrides InitialSlope
%
%Example usage:
%
% 	DAE = BJTdiffpairSchmittTrigger();
%	qss = QSS(DAE); qss = feval(qss.solve, qss); sol = feval(qss.getsolution, qss);
%	trans = ODEXY(DAE);
%	trans = feval(trans.solve, trans, sol, 0, 1e-5, 10e-3);
%	feval(trans.plot, trans); % bad results - takes timesteps that are too large
%	tranparms.MaxStep = 1e-5;
%	trans = ODEXY(DAE, [], tranparms);
%	trans = feval(trans.solve, trans, sol, 0, 1e-5, 10e-3);
%	feval(trans.plot, trans); % good results
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Type "help MAPPlicense" at the MATLAB/Octave prompt to see the license      %
%% for this software.                                                          %
%% Author: J. Roychowdhury, 2009/sometime                                      %
%% Copyright (C) 2008-2020 Jaijeet Roychowdhury <jr@berkeley.edu>. All rights  %
%%               reserved.                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%
	if (nargin > 3) || (nargin < 1)
		fprintf(2,'ODEXY: error: too many or too few arguments.\n');
		help('ODEXY');
		return;
	end
	ODEXYobj = transient_skeleton(DAE); % set up basic structure and useful
					  % functions

	% usage and name strings
	ODEXYobj.Usage = help('ODEXY'); 
	% name is a function

	%
	ODEXYobj.name = @name;
	ODEXYobj.solve = @ODEXYsolve;

	if (nargin < 2) || isempty(TRmethod)
		ODEXYobj.TRmethod = 'ode15s';
	else
		ODEXYobj.TRmethod = TRmethod;
	end

	% Check and set the solver method
	switch ODEXYobj.TRmethod;
	  case {'ode15s'}
	     ODEXYobj.solvefunc = @ode15s;
	  case {'ode45', 'ode23'}
	     fprintf(2, 'ODEXY: using solver %s (works ONLY for ODEs: ie, q(x) MUST EQUAL -x)\n', ODEXYobj.TRmethod);
	     x = rand(feval(DAE.nunks, DAE), 1);
	     q = feval(DAE.q, x, DAE);
	     if norm(x+q) ~= 0
	     	fprintf(2, 'ODEXY: ERROR: q(x) != -x for x=\n'); x
	     	fprintf(2, 'Please fix your DAE or change your transient solver method (to, eg, ode15s)\n');
		ODEXYobj = []; return;
	     end
	     if strcmp(ODEXYobj.TRmethod, 'ode45')
	     	ODEXYobj.solvefunc = @ode45;
	     elseif strcmp(ODEXYobj.TRmethod, 'ode23')
	     	ODEXYobj.solvefunc = @ode23;
	     end

	  otherwise
	     fprintf(2, 'ODEXY: ERROR: unsupported solver method %s\n', ODEXYobj.TRmethod);
		ODEXYobj = []; return;
	end


	if (nargin >= 3)
		ODEXYobj.tranparms = tranparms;
	end

	% Set up ode options
	ODEXYobj.ODEoptions = odeset('RelTol', 1.0e-4, 'Abstol', 1.0e-9);
	if isfield(ODEXYobj.tranparms, 'RelTol')
		ODEXYobj.ODEoptions = odeset(ODEXYobj.ODEoptions, ...
			'RelTol', ODEXYobj.tranparms.RelTol);
	end
	if isfield(ODEXYobj.tranparms, 'AbsTol')
		ODEXYobj.ODEoptions = odeset(ODEXYobj.ODEoptions, ...
			'AbsTol', ODEXYobj.tranparms.AbsTol);
	end
	if isfield(ODEXYobj.tranparms, 'MaxStep')
		ODEXYobj.ODEoptions = odeset(ODEXYobj.ODEoptions, ...
			'MaxStep', ODEXYobj.tranparms.MaxStep);
	end
	if isfield(ODEXYobj.tranparms, 'InitialStep')
		ODEXYobj.ODEoptions = odeset(ODEXYobj.ODEoptions, ...
			'InitialStep', ODEXYobj.tranparms.InitialStep);
	end
	if isfield(ODEXYobj.tranparms, 'InitialSlope')
		ODEXYobj.ODEoptions = odeset(ODEXYobj.ODEoptions, ...
			'InitialSlope', ODEXYobj.tranparms.InitialSlope);
	end
	if isfield(ODEXYobj.tranparms, 'MaxOrder')
		ODEXYobj.ODEoptions = odeset(ODEXYobj.ODEoptions, ...
			'MaxOrder', ODEXYobj.tranparms.MaxOrder);
	end
	if isfield(ODEXYobj.tranparms, 'BDF')
		ODEXYobj.ODEoptions = odeset(ODEXYobj.ODEoptions, ...
			'BDF', ODEXYobj.tranparms.BDF);
	end
	ODEXYobj.BEforFirstStep = 1;
	if isfield(ODEXYobj.tranparms, 'BEforFirstStep')
		if 0 == ODEXYobj.tranparms.BEforFirstStep
			ODEXYobj.BEforFirstStep = 0;
		end
	end

	% Set up functions for the solver

	% f(t,y) 
	if 1 == ODEXYobj.f_takes_u
		ODEXYobj.f = @(t, y) ODEXYobj.DAE.f(y, ...
			ODEXYobj.DAE.utransient(t, ODEXYobj.DAE), ODEXYobj.DAE);
		% Jacobian of f: dfdy(t,y) 
		ODEXYobj.dfdy = @(t, y) ODEXYobj.DAE.df_dx(y, ...
			ODEXYobj.DAE.utransient(t, ODEXYobj.DAE), ODEXYobj.DAE);
	else
		ODEXYobj.f = @(t, y) ODEXYobj.DAE.f(y, ODEXYobj.DAE) + ...
			ODEXYobj.B*ODEXYobj.DAE.utransient(t, ODEXYobj.DAE);
		% Jacobian of f: dfdy(t,y) 
		ODEXYobj.dfdy = @(t, y) ODEXYobj.DAE.df_dx(y, ODEXYobj.DAE);
	end
	% "Mass matrix": Jacobian of -q()
	ODEXYobj.Mass = @(t, y) -ODEXYobj.DAE.dq_dx(y, ODEXYobj.DAE);

	ODEXYobj.ODEoptions = odeset(ODEXYobj.ODEoptions, ...
				     'Mass', ODEXYobj.Mass, ...
				     'Jacobian', ODEXYobj.dfdy);
end
% end of ODEXY "constructor"

function out = name(ODEXYobj)
	out = sprintf('%s ODEXY solver', ODEXYobj.TRmethod);
end
% end name

function ODEXYobjOUT = ODEXYsolve(ODEXYobj, xinitcond, tstart, tstep, tstop, store_Jacobians)
%function ODEXYobjOUT = ODEXYsolve(ODEXYobj, xinitcond, tstart, tstep, tstop, store_Jacobians)

	if nargin < 6
		store_Jacobians = 0;
	end

	if 1 == ODEXYobj.BEforFirstStep % run BE for 1 small step
		BEobj = LMS(ODEXYobj.DAE);
		BEobj = feval(BEobj.solve, BEobj, xinitcond, tstart, tstep/10, tstep/10, store_Jacobians);
		[ODEXYobj.tpts, ODEXYobj.vals] = feval(BEobj.getsolution, BEobj);
		tstart = ODEXYobj.tpts(BEobj.timeptidx);
		xinitcond = ODEXYobj.vals(:, BEobj.timeptidx);
		xdotinitcond = (ODEXYobj.vals(:, BEobj.timeptidx) - ODEXYobj.vals(:, 1)) ...
				/ ( ODEXYobj.tpts(BEobj.timeptidx) - ODEXYobj.tpts(1) );
		ODEXYobj.ODEoptions = odeset(ODEXYobj.ODEoptions, 'InitialSlope', xdotinitcond);
	end

	[ts, ys] = feval(ODEXYobj.solvefunc, ODEXYobj.f, [tstart tstop], ...
			xinitcond, ODEXYobj.ODEoptions);
	ODEXYobj.tpts = [ODEXYobj.tpts(1:(end-1)), ts.'];
	ODEXYobj.vals = [ODEXYobj.vals(:,1:(end-1)), ys.'];

	ODEXYobj.timeptidx = length(ODEXYobj.tpts);

	if 1 == store_Jacobians
		for i = 1:length(ODEXYobj.tpts)
			t = ODEXYobj.tpts(i);
			x = ODEXYobj.vals(:,i);
			u = feval(ODEXYobj.DAE.utransient, t, ODEXYobj.DAE);
			ODEXYobj.Cs{i} = feval(ODEXYobj.DAE.dq_dx, x, ODEXYobj.DAE);
			ODEXYobj.Gs{i} = feval(ODEXYobj.DAE.df_dx, x, u, ODEXYobj.DAE);
			ODEXYobj.Gus{i} = feval(ODEXYobj.DAE.df_du, x, u, ODEXYobj.DAE);
		end
	end

	ODEXYobjOUT = ODEXYobj;
	ODEXYobjOUT.solvalid = 1;
end
% end of solve
