function SHOOTparms = defaultShootingParms()
	% shooting NR parameters
	%
	SHOOTparms.periodicity_reltol = 1e-6; % not used at the moment
	SHOOTparms.periodicity_abstol = 1e-12; % not used at the moment
	%
	SHOOTparms.Nsteps = 100; % try 100 time-steps per period
	%
	SHOOTparms.dbglvl = 2; 
	% dbglvl (default 0; can be -1, 0, 1 and 2)
	%		-1: print nothing (not even errors)
	%		0: only print errors.
	%		1: print . per shooting NR iteration to indicate progress
	%			- suppress all transient output except for errors
	%		2: print transient progress output as well
	%			- at transient dbglvl 1
	TRmethods = LMSmethods();
	%SHOOTparms.TRmethod = TRmethods.BE;
	SHOOTparms.TRmethod = TRmethods.GEAR2;

	%
	SHOOTparms.tranparms = defaultTranParms();
	SHOOTparms.tranparms.trandbglvl = SHOOTparms.dbglvl-1;

	%
	NRparms = SHOOTparms.tranparms.NRparms; % not defaultNRparms
		% shooting tolerances should be looser
		% than the tolerances of the transient underlying it.
	NRparms.reltol = 5*NRparms.reltol;
	NRparms.abstol = 5*NRparms.abstol;
	NRparms.residualtol = 5*NRparms.residualtol;
	NRparms.limiting = 0; % transient does its own limiting
	SHOOTparms.NRparms = NRparms;
