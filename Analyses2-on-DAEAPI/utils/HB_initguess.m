function Xinitguess_Nn = HB_initguess(DAE, vecf0, vecM, DCguess, transient_initcond, DCorTransient, cycles, tsteps_per_cycle, cycles_to_skip, harmonics_to_keep, doplot)
%function Xinitguess_Nn = HB_initguess(DAE, vecf0, vecM, DCguess, ...
%                                      transient_initcond, DCorTransient, ...
%                                      cycles, tsteps_per_cycle, ...
%                                      cycles_to_skip, harmonics_to_keep, ...
%                                      doplot)
%Returns an initial guess for the HB Newton-Raphson by running a transient 
%	for a few cycles and using its Fourier components for the guess.
%
%Arguments:
%
%- DAE: the DAE. It must have uHB set up before this function is called.
%- vecf0: Frequencies at which HB will be run.
%- vecM: number of positive harmonics for each fundamental frequency.
%
%---------- the remaining arguments are optional:
%- DCguess: NR guess for DC analysis, if conducted. If not supplied (or [])
%- transient_initcond: initial condition for the transient analysis. If not
%	supplied or [], defaults to the DC solution (calculated by this routine).
%- DCorTransient: a string. If 'DC', then only the DC operating point will
%   be used for the HB initial condition (ie, no transient will be run),
%	and vecf0, vecM, transient_initcond, etc. will be ignored.  The default is to
%	use a transient.
%- cycles: number of 1/f0 periods for which the transient will be run.
%	Must be >= 2.  Using an even number of cycles is highly recommended.
%	Defaults to 4 if not specified or [].
%- tsteps_per_cycle: number of timesteps per cycle for the transient analysis.
%	Defaults to 30 if not specified or [].
%- cycles_to_skip: how many cycles from the start to skip for the fft.
%   Defaults to 0 if not specified or [].
%- harmonics_to_keep: how many harmonics of vecf0 to use from the transient.
%	Should be <= vecM. Defaults to M if not specified or [].
%- doplot: if 1, plot transient results (all outputs). Default is not to plot.
%	
%Changelog:
%2015/01/07, JR:
%   - added the DCguess argument, to be passed to QSS.solve if specified.
%   - was using tstop instead of unif_tpts(end), leading to 
%     junk being returned when cycles_to_skip > 0, leading to
%     test_LCtanhOsc_HB always converging to the DC operating point.
%   - Also, there was complex conjugacy error when harmonics_to_keep
%     was less than M. Both fixed, confirmed test_LCtanhOsc_HB 
%     now converging properly.
%
%Changelog:
%2017/04/06, AG:
%   - We would like to use the same script to generate an initial guess for 
%     multi-tone Harmonic balance. However, it seems to have been hardcoded
%     for 1 fundamental frequency.
%   - SETTING UP TRANSIENT INPUTS:
%     [TODO]
%            
%Author: Jaijeet Roychowdhury <jr@berkeley.edu> 2012/06/08
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J. Roychowdhury.
% Copyright (C) 2008-2012 Jaijeet Roychowdhury <jr@berkeley.edu>. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

	if nargin < 4 || isempty(DCguess) 
		DCguess = [];
	end

	if nargin < 5 || isempty(transient_initcond) 
		transient_initcond = [];
	end

	if nargin < 6 || isempty(DCorTransient)
		DCorTransient = 'transient';
	end

	if nargin < 7 || isempty(cycles)
		cycles = 4;
	end

	if nargin < 8 || isempty(tsteps_per_cycle)
		tsteps_per_cycle = 30;
	end

	if nargin < 9 || isempty(cycles_to_skip)
		cycles_to_skip = 0;
	end

	if nargin < 10 || isempty(harmonics_to_keep)
		harmonics_to_keep = vecM;
	end

	if nargin < 11 || isempty(doplot)
		doplot = 0;
	end


	% get the input harmonics for HB
	nu = length(DAE.uHBfunc_updates);
	if nu > 0
		U_twoD = feval(DAE.uHB, vecf0, vecM, DAE);
	else
		U_twoD = [];
	end

	% run a DC op: if DC only requested, or if transient_initial cond not specified
	if ( (0 == length(transient_initcond)) || (1==strcmp(DCorTransient,'DC')) )
		fprintf(1,'running an initial DC...\n');
		if ni > 0
			DAE = feval(DAE.set_uQSS, real(U_twoD(:,1)), DAE);
		end
		qss = QSS(DAE);
        if isempty(DCguess)
		    qss = feval(qss.solve, qss);
        else
		    qss = feval(qss.solve, DCguess, qss);
        end
		DCsol = feval(qss.getsolution, qss);
		fprintf(1,'\n...DC solution obtained.\n');
	end

	% set up N for HB
	N=2*vecM+1; 

	if (1 ~= strcmp(DCorTransient,'DC')) % use a short transient to set up HB initial guess
		% set up transient inputs corresponding to the HB input
		utargs.U_twoD = U_twoD;
		utargs.f = vecf0;
		DAE = feval(DAE.set_utransient, @utfunc_from_U, utargs, DAE);

		% set up a transient for a few cycles
		T = 1./vecf0;
		if (0 == length(transient_initcond))
			xinit = DCsol;
		else
			xinit = transient_initcond;
		end
		tstart = 0;
		tstep  = min(T)/tsteps_per_cycle;   % Need to consider the largest frequency (smallest period) for the
                                                    % transient time-step
		tstop  = max(T)*cycles; % Similarly, the smallest frequency has to be considered for figuring out the
                                        % total length of the transient simulation... Getting enough of all the
                                        % frequencies

		fprintf(1,'\nrunning a transient of %d cycles using GEAR2 to obtain an initial guess for HB...\n', cycles);
		% run the transient
		tran = run_transient_GEAR2(DAE, xinit, tstart, tstep, tstop);
		[tpts, transol] = feval(tran.getsolution, tran);
		if 0 ~= doplot
			souts = StateOutputs(DAE);
			feval(tran.plot, tran, souts);
		end

		firstidx = min(find(tpts >= T*cycles_to_skip)); % index of first tpt >= T*cycles_to_skip;
		new_tpts = tpts(firstidx:end) - tpts(firstidx);
		new_transol = transol(:,firstidx:end);


		% interpolate the transient solution on to a uniform grid
		N_ts = 2*(cycles-cycles_to_skip)*N; % compute enough harmonics to cover vecM*vecf0
		unif_tpts = (0:(N_ts-1))/(N_ts-1)*(tpts(end)-tpts(firstidx));
		interpolated_vals = interp1(new_tpts, new_transol.', unif_tpts).';
        % this plot for debugging
        %for unkidx = 1:size(transol,1)
        %    plot(unif_tpts+tpts(firstidx), interpolated_vals(unk_idx, :), 'o-');
        %end

		% run an fft on it to find its fourier components approximately
                % [TODO]: This probably needs to be changed
		fcoeffs = 1/N_ts*fft(interpolated_vals.').';
		%fcoeffs = N_ts*fft(interpolated_vals.').';

		% find the fourier component corresponding most closely to freq vecf0
		  % how to figure out the fundamental harmonic for the above fcoeffs:
		  % given a set of fourier coeffs X with fundamental freq f,
		  % x = N*ifft(X) produces N equally space time-domain samples over 
		  % 1 period = 1/f; ie, the spacing between the time-domain samples
		  % is 1/(N*f). Hence we have unif_tpts(1) = 1/(N*f), or
		  %	f = 1/(N*unif_tpts(1)) = 1/tstop.

		f_fundamental = 1/unif_tpts(end);
		dc_component = fcoeffs(:,1);
		n = size(dc_component,1);
		Xinitguess_twoD(:,1) = dc_component;
		
		% find approximate value for each harmonic of vecf0
		for i=1:vecM
			idx = round(1 + i*vecf0/f_fundamental);
			if idx < 2
				fprintf(2,'error in transient harmonic component estimation: idx=%d<2 for harmonic i=%d\n', idx,i);
				return;
			end
			if i <= harmonics_to_keep 
				Xinitguess_twoD(:,i+1) = fcoeffs(:,idx);
			    Xinitguess_twoD(:,N-i+1) = conj(fcoeffs(:,idx));
			else
				Xinitguess_twoD(:,i+1) = zeros(n,1);
			    Xinitguess_twoD(:,N-i+1) = zeros(n,1);
			end
		end

		% set up Xinitguess_Nn 
		Xinitguess_Nn = twoD_to_Nn(Xinitguess_twoD);
		fprintf(1,'\n...transient-based HB initial guess obtained.\n');
	else % otherwise, use the DC solution as HB initial guess
		% set up the initial guess for HB as the DC solution
		n = length(DCsol);
		DCsol_twoD = sparse(n, N); DCsol_twoD(:,1) = DCsol;
		Xinitguess_Nn = full(twoD_to_Nn(DCsol_twoD));
		fprintf(1,'\nUsing the DC solution as the initial guess for HB.\n');
	end
end
