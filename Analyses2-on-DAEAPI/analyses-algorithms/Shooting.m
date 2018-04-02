function SHOOTobj = Shooting(DAE, SHOOTparms, isOsc) % DAE=DAEAPIv6.2
%function SHOOTobj = Shooting(DAE, SHOOTparms, isOsc) % DAE=DAEAPIv6.2
%Author: Jaijeet Roychowdhury <jr@berkeley.edu>, 2011/12/14-20
%
%SHOOTparms.TRanparms - passed to SHOOT
%SHOOTparms.TRmethod - passed to SHOOT
%SHOOTparms.T - time period (this could be arg to Solve?)
%SHOOTparms.periodicity_reltol - reltol for checking periodicity
%SHOOTparms.periodicity_abstol - abstol for checking periodicity
%note: DAE should be set up with a T-periodic utransient 
%	- should check it here, using some random technique
%
%Shooting: trannsient simulation object
%	- calls SHOOT, SHOOTforTabulatedLTV and NR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to define the state transition function Phi(t;t0,x0) and its derivative.
% the derivative should probably be implemented as an optional part of SHOOT
% (this should also be useful for transient sensitivities).
% - done: SHOOTforTabulatedLTV
%	- see also ~/.../xournal-notes/2011-12-13-Note-22-32_LTVDAE.xoj
%
% next step: write down the theory in detail, especially with regard to the
% propagation of dPhi/dx0 (and maybe dPhi/dt, for oscillators).
% 
% - done: ~/.../xournal-notes/2011-12-12-Note-17-39_shooting.xoj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Example usage: TODO
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Type "help MAPPlicense" at the MATLAB/Octave prompt to see the license      %
%% for this software.                                                          %
%% Author: Jaijeet Roychowdhury <jr@berkeley.edu>, 2011/12/14-20               %
%% Copyright (C) 2008-2020 Jaijeet Roychowdhury <jr@berkeley.edu>. All rights  %
%%               reserved.                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
	if (nargin > 3) || (nargin < 1)
		fprintf(2,'Shooting: error: too many or too few arguments.\n');
		help('Shooting');
		return;
	end

	% usage and name strings
	SHOOTobj.Usage = help('Shooting'); 
	% name is a function

	% data/precomputation
	SHOOTobj.version = 'DAEAPIv6.2';
	SHOOTobj.DAEversion = 'DAEAPIv6.2';
	%
	%
	SHOOTobj.DAE = DAE;
	SHOOTobj.solvalid = 0;
	SHOOTobj.tpts = [];
	SHOOTobj.vals = [];
	SHOOTobj.solution = [];

	if nargin < 2 || isempty(SHOOTparms)
		SHOOTparms = defaultShootingParms();
	end
	SHOOTobj.parms = SHOOTparms;

	if nargin < 3
		isOsc = 0;
	end
    
    if 1 == isOsc
        SHOOTparms.NRparms.MPPINR_use_pinv = 1; 
    end

	SHOOTobj.isOsc = isOsc;


	%FIXME: check that DAE's utransient is defined and periodic (but we don't have T yet)

	SHOOTobj.tran = LMS(DAE, SHOOTparms.TRmethod, SHOOTparms.tranparms);

	% function handles
	SHOOTobj.solve = @SHOOTsolve;
	SHOOTobj.getsolution = @SHOOTgetsolution;
	SHOOTobj.plot = @SHOOTplot;
	SHOOTobj.updateDAE = @SHOOTupdateDAE;
	SHOOTobj.STF_dSTF = @STF_dSTF;
	if SHOOTparms.dbglvl > 2
		SHOOTobj.STF_dSTF = @STF_dSTF;
		SHOOTobj.gJfunc_for_NR = @gJfunc_for_NR;
	end
end
% end of SHOOT "constructor"

function SHOOTout = SHOOTupdateDAE(DAE, SHOOTobj)
	SHOOTobj.DAE = DAE;
	SHOOTobj.tran = feval(SHOOTobj.tran.updateDAE, DAE, SHOOTobj.tran);
	SHOOTout = SHOOTobj;
end
% end of SHOOTupdateDAE

function SHOOTout = SHOOTsolve(SHOOTobj, xinit, T, tstart)
	% T is a guess if isOsc == 1
	SHOOTobj.T = T;
	if nargin < 4
		tstart = 0;
	end

	% somewhat kludgy way of passing info to gJfunc_for_NR
	SHOOTobj.tinit = tstart;

	if 1 == SHOOTobj.isOsc
		Tr = 1; % relative T
		xinit(end+1,1) = Tr; % the unknown is Tr, not T. T = Tr*SHOOTobj.T
	end

	% run Newton-Raphson
	if SHOOTobj.parms.dbglvl > 1 
		fprintf(2,'[start shooting NR...\n');
	end
	[solution, iters, success] = NR(@gJfunc_for_NR, [], xinit, SHOOTobj, SHOOTobj.parms.NRparms);
	if SHOOTobj.parms.dbglvl > 1 
		fprintf(2,'\n...shooting NR completed]\n');
	end
	if 1 ~= success
		if SHOOTobj.dbglvl > -1 
			fprintf(2,'\nShooting: NR failed\n');
		end
		SHOOTobj.solvalid = 0; SHOOTobj.tpts = []; SHOOTobj.vals = []; SHOOTobj.solution = [];
		SHOOTout = SHOOTobj;
		return;
	end

	if SHOOTobj.parms.dbglvl > 0 
		fprintf(2,'\nShooting: NR succeeded in %d iterations.\n', iters);
	end

	if 1 == SHOOTobj.isOsc 
		Tr = solution(end);
		T = Tr*SHOOTobj.T; % last entry of solution is Tr
		solution = solution(1:(end-1));
	end

	% run a final transient simulation to get the full track
	tend = tstart + T; 
	tstep = T/SHOOTobj.parms.Nsteps;
	if SHOOTobj.parms.dbglvl > 1 
		fprintf(2,'\nstart final transient run[\n');
	end
	[xend, dxend_dx0, success, tpts, vals, ltvtran] = STF_dSTF(solution, tstart, tstep, tend, SHOOTobj);
	if SHOOTobj.parms.dbglvl > 1 
		fprintf(2,'final transient completed]\n');
		normsol = norm(xend);
		%normerr = norm(xend-vals(:,1));
		normerr = norm(xend-solution);
		fprintf(2,'periodicity error (x(0) vs x(T)): absolute: %g, relative: %g\n', ...
			normerr, normerr/(normsol+1e-18));
	end
	if 1 ~= success
		error('Shooting: weird error: shooting NR succeeded, but final transient didn''t.');
	end
	if 1 == SHOOTobj.isOsc
		SHOOTobj.Tinitial = SHOOTobj.T;
		SHOOTobj.Tr = Tr;
	end
	SHOOTobj.T = T;
	SHOOTobj.solution = solution;
	SHOOTobj.tpts = tpts;
	SHOOTobj.vals = vals;
	SHOOTobj.MonodromyMatrixAtSolution = dxend_dx0; % if isOsc==1, this is the augmented monodromy matrix: last col is dx_end/dT
	SHOOTobj.LTVtran = ltvtran; % LMSLTVtran object with G(t), C(t), etc. set up at solution
	SHOOTobj.solvalid = 1;
	SHOOTout = SHOOTobj;
end
% end of SHOOTsolve

function [g, J, success] = gJfunc_for_NR(x, SHOOTobj)
	tstart = SHOOTobj.tinit; % set up in solve
	T = SHOOTobj.T;
	if 1 == SHOOTobj.isOsc
		Tr = x(end);
		T = SHOOTobj.T*Tr;
		x = x(1:(end-1));
	end
	tend = tstart + T; 
	tstep = T/SHOOTobj.parms.Nsteps;
	[xend, dxend_dx0, success] = STF_dSTF(x, tstart, tstep, tend, SHOOTobj);

	% shooting: g(x) = Phi(t0+T; t0, x0) - x0;
	g = xend - x;
	J = dxend_dx0 - eye(size(dxend_dx0));
	% the above works for isOsc == 1 as well, since
	% eye(n,n+1) produces zeros for the rightmost column
	% ie, no phase condition necessary.
end
% end of gJfunc_for_NR

function [xend, dxend_dx0, success, tpts, vals, ltvtran] = STF_dSTF(xinit, tstart, tstep, tend, SHOOTobj)
	tran = SHOOTobj.tran;
	store_Jacobians = 1;
	% run transient, store Jacobians
	if SHOOTobj.parms.dbglvl > 1 
		fprintf(2,' <start transient+dx/dx0 run...\n\t');
	end
	if tend <= tstart ||tstart + tstep >= tend 
		fprintf(2,' error: tend <= tstart or tstart + tstep >= tend ...\n\t');
	end
	tran = feval(tran.solve, tran, xinit, tstart, tstep, tend, store_Jacobians);
	if 1 ~= tran.solvalid
		success = tran.solvalid;
		xend = []; dxend_dx0 = []; tpts = []; vals = [];
		% FIXME: throw an error message
		return;
	end
	[tpts, vals, jacobians] = feval(tran.getsolution, tran);

	xend = vals(:, end);
	if 1 == SHOOTobj.isOsc
		p = tran.TRmethod.order;
		last_pp1_ts(1) = tpts(end);
		last_pp1_ts(2:(p+1)) = tpts(end-1:-1:end-p);
		last_pp1_xs(:,1) = vals(:,end);
		last_pp1_xs(:,2:(p+1)) = vals(:,end-1:-1:end-p);
	end

	jacobians.tpts = tpts;
	jacobians.Sus = []; % no input to LTV system, just initial cond
	nunks = feval(SHOOTobj.DAE.nunks, SHOOTobj.DAE);
	initdxdx0 = eye(nunks, nunks);

	% run LTV transient to get dx/dx0
	ltvtran = LMSforTabulatedLTV(jacobians, SHOOTobj.parms.TRmethod, ...
		SHOOTobj.parms.tranparms);
	ltvtran = feval(ltvtran.solve, ltvtran, initdxdx0);
	if SHOOTobj.parms.dbglvl > 1 
		fprintf(2,'  ...transient+dx/dx0 run completed>\n');
	end
	if 1 ~= ltvtran.solvalid
		success = ltvtran.solvalid;
		xend = []; dxend_dx0 = []; tpts = []; vals = [];
		% FIXME: throw an error message
		return;
	end
	[tpts, vals2] = feval(ltvtran.getsolution, ltvtran);
	dxend_dx0 = vals2{end};
	if 1 == SHOOTobj.isOsc
		% we want d/dT x(T)
		% d/dT x(T) = x(T+delta)-x(T)/delta = d/dt x(t) @ t=T
		%
		% what we really want is d/dTr x(T), 
		%	with Tr = T/SHOOTobj.T => dT = SHOOTobj.T * dTr
		% 	d/dTr x(T) = SHOOTobj.T * d/dt x(t) @ t=T
		lastcol = SHOOTobj.T * feval(tran.TRmethod.ddtApproxATtn, last_pp1_ts, last_pp1_xs);
		dxend_dx0 = [dxend_dx0, lastcol];
	end
	success = 1;
end
% end of STF_dSTF

%function [solution, tpts, vals, T, MonodromyMatrixAtSolution, LTVtranObjAtSolution] = SHOOTgetsolution(SHOOTobj)
function shootsol = SHOOTgetsolution(SHOOTobj)
% returns solution over one period, including first and last points
	if SHOOTobj.solvalid < 1
		fprintf(2,'Shooting: getsolution: run solve successfully first!\n');
		shootsol = [];
	else
		shootsol.solution = SHOOTobj.solution;
		shootsol.tpts = SHOOTobj.tpts;
		shootsol.vals = SHOOTobj.vals;
		shootsol.T = SHOOTobj.T; % useful mainly for isOsc == 1.
		shootsol.MonodromyMatrixAtSolution = SHOOTobj.MonodromyMatrixAtSolution;
		shootsol.LTVtranObjAtSolution = SHOOTobj.LTVtran;
	end
end
% end SHOOTgetsolution

function [figh, onames, colindex] = SHOOTplot(SHOOTobj, varargin)
	DAE = SHOOTobj.DAE; tpts = SHOOTobj.tpts; vals = SHOOTobj.vals;
	time_units = DAE.time_units;
	[figh, onames, colindex] = transientPlot(DAE, tpts, vals, time_units, varargin{:});
	if 1 == SHOOTobj.isOsc
		ttlstr = sprintf('%s: periodic self-oscillatory solution via shooting using %s\nT=%gs, f=%gHz', ...
			feval(SHOOTobj.DAE.daename,SHOOTobj.DAE),  ...
			SHOOTobj.parms.TRmethod.name, SHOOTobj.T,1/SHOOTobj.T);
	else
		ttlstr = sprintf('%s: periodic solution via shooting using %s', ...
			feval(SHOOTobj.DAE.daename,SHOOTobj.DAE), ...
			SHOOTobj.parms.TRmethod.name);
	end
	title(ttlstr);
end
% end of function SHOOTplot
