function HBobj = HB(DAE, isOsc, NRparms_for_HB)
%function HBobj = HB(DAE, isOsc=0, NRparms_for_HB=[default_NR_parms])
%1-tone harmonic balance analysis.
%
%Arguments:
%- DAE: the DAE to run HB on. It should have HB inputs defined for its inputs.
%- isOsc is optional. It should be set to 1 if the system is autonomous. 
%   Default if not specified is 0 == non-autonomous
%- NRparms is optional. These are the NR parameters to be used for solving
%   the HB equations using Newton-Raphson. If specified, it should
%   be a structure with the fields specified in utils/defaultNRparms.m.
%   If not specified, defaultNRparms are used.
%
%Methods (function handles) of HBobj:
%
%- HBobj = solve(HBobj, Xinitguess_Nn, M, f0)
%  - Xinitguess_Nn: the initial guess for NR, in Nn format
%    - if an oscillator, this should NOT include the frequency
%      guess as an unknown - that should be supplied as the f0
%      argument.
%  - M: the number of 1-sided harmonics (excluding DC) to
%       be used in the 1-tone HB analysis.
%  - f0: the fundamental frequency for non-autonomous HB analysis.
%      For autonomous systems (ie, isOsc == 1), f0 is used
%      as an initial guess for NR.
%
%- hbsolStruct = getsolution(HBobj)
%  - returns the HB solution (after successful solve) in a structure
%    hbsolStruct with the following fields:
%    .frequency: fundamental frequency of the HB analysis. If isOsc==1, this
%                was solved for as part of HB.
%    .X_twoD: 2-D (n x N matrix) form of the HB solution - the Fourier
%                coefficients in standard fft order.
%    .x_twoD: 2-D (n x N matrix) form of the HB solution, converted to the time
%                domain (equispaced samples over the period).
%    .tpts: the timepoints corresponding to the columns of x_twoD
%    .U_twoD: 2-D (n x N matrix) form of the inputs to the DAE system - the
%                Fourier coefficients in standard fft order.
%    .u_twoD: 2-D (n x N matrix) form of the inputs to the DAE system -
%                time-domain samples.
%    .Jacobian_Nne_Nn: the Jacobian of the HB system H(X) (corresponding to Nn
%                format for both H and X). This is of size 
%                - N*neqns x N*nunks if isOsc == 0,
%                - N*neqns x (N*nunks+1) if isOsc == 1.
%
%- [figh, onames, colindex] = plot(HBobj, [HBobj, stateoutputs, ...
%                          lgndprefix, linetype, figh, legends, colindex])
%  - all arguments except the first are optional
%  - if additional arguments are not specified, bar-plots all the harmonics of
%    the DAE-defined outputs.
%  - if stateoutputs is specified (as a StateOutputs object, see help
%    StateOuputs), all harmonics of the specified outputs are plotted.
%  - to understand usage with the additional arguments and the outputs, see
%    help plot_FD_TD and help transientplot.
%
%- HBobj = compute_sensitivities_direct(pObj, HBobj)
%  - computes sensitivities using a direct (ie, not adjoint) algorithm,
%    producing a potentially large, dense sensitivity matrix S_HB, which
%    is stored as HBobj.S_HB. solve() should have been called successfully
%    first for the results to be meaningful.
%
%- [TODO] HBobj = compute_sensitivities_adjoint(pObj, HB_outputs_obj, HBobj) -
%         yet to be implemented
%
%- S_HB = getSensitivities(HBobj)
%  - returns the sensitivity matrix S_HB set up by 
%    compute_sensitivities_direct(). Adjoint support will be added later.
%
%- [TODO] plotSensitivities(HBobj) - not implemented yet.
%
%TODO: support oversampling
%
%TODO: usage example
%

%% Author: J. Roychowdhury, 2012/{05/25-06/05}, 2017/03/19-22                  %

%% Changelog:
%   2017/03/19: added direct_sensitivities - JR.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Type "help MAPPlicense" at the MATLAB/Octave prompt to see the license      %
%% for this software.                                                          %
%% Copyright (C) 2008-2020 Jaijeet Roychowdhury <jr@berkeley.edu>. All rights  %
%%               reserved.                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % store arguments
    HBobj.DAE = DAE;
    if nargin > 1
        HBobj.isOsc = isOsc;
    else
        HBobj.isOsc = 0;
    end
    if nargin > 2
        HBobj.parms.NRparms = NRparms_for_HB;
    else
        HBobj.parms.NRparms = defaultNRparms();
    end

    if 1 == HBobj.isOsc
        HBobj.parms.NRparms.MPPINR_use_pinv = 1; % use pinv() based solution
    end

    % set up solution storage
    HBobj.solution = [];
    HBobj.Jacobian_at_sol = [];
    HBobj.solvalid = 0;

    % set up function pointers for the analysis
    HBobj.solve = @solve_HB;
    HBobj.plot = @plot_HB;
    HBobj.getsolution = @getsolution_HB;

    % parametric sensitivity analysis functions
    HBobj.compute_sensitivities_direct = @compute_sensitivities_direct;
    HBobj.getSensitivities = @getSensitivities_HB;

end % HB constructor

function [g, J, evalsuccess] = H_HB(X_Nn, HBobj)
%function [g, J, evalsuccess] = H_HB(X_Nn, HBobj)
%This function (used directly by NR.m) returns
%the H(X_Nn) function for HB, and its Jacobian.
%See the HB documentation in doc/.
    % global HBdebuglvl; % REMOVE
    f = HBobj.f0; % freq. known a priori
    if 1 == HBobj.isOsc
        fr = X_Nn(end,1); % last entry of the unknown vector is the relative frequency fr

        % hack to keep numerical errors in imag components from going out of hand
        if 1 == 0
            % DOES NOT HELP: ENABLING THIS DOES NOT HELP NR
            fr = real(fr);
        end
        % end hack to keep numerical errors in imag components from going out of hand

        f = HBobj.f0*fr; % PUT THE SEMICOLON BACK

        X_Nn = X_Nn(1:(end-1),1);

        % re-evaluate U at this frequency
        HBobj.U_Nni = setup_UNni(HBobj.DAE, f, HBobj.M, HBobj.N);
    end

    % hack to clean up complex conjugacy errors in X_Nn that may have been introduced somehow by NR
    if 1 == 0
        % ENABLING THIS TOTALLY BREAKS NR.
        n = feval(HBobj.DAE.nunks, HBobj.DAE);
        N = HBobj.N; M = HBobj.M;
        X_twoD = Nn_to_twoD(X_Nn, n, N);
        X_twoD(:,end:-1:(end-M+1)) = conj(X_twoD(:,2:(M+1)));
        % DC component
        X_twoD(:,1) = real(X_twoD(:,1));
        X_Nn = twoD_to_Nn(X_twoD);
    end
    % end hack to clean up complex conjugacy errors in X_Nn that may have been introduced somehow by NR


    %{
    if 1 == HBdebuglvl
        n = feval(HBobj.DAE.nunks, HBobj.DAE)
        N = HBobj.N;
        Nn_to_twoD(X_Nn, n, N)
    end
    %}


    % these two functions (in utils) are where the real work is done. See
    % comments within the files, and also the HB documentation. They are kept
    % separately because there are standalone tests to verify their
    % correctness.
    [F_Nne, Q_Nne, dFdX, dQdX] = FQJ_HB(X_Nn, HBobj.U_Nni, HBobj);
    %{
    if 1 == HBdebuglvl
        ne = feval(HBobj.DAE.neqns, HBobj.DAE)
        F_Nne_twoD = Nn_to_twoD(F_Nne, ne, N)
        Q_Nne_twoD = Nn_to_twoD(Q_Nne, ne, N)
    end
    %}
    [OmegaQ_Nne, dOmegaQdX] = OmegaQ_HB(Q_Nne, dQdX, f/HBobj.f0, HBobj); % isOsc==1 case is handled within OmegaQ_HB
    %{
    if 1 == HBdebuglvl
        OmegaQ_Nne_twoD = Nn_to_twoD(OmegaQ_Nne, ne, N)
    end
    %}

    % g = H_HB = OmegaQ_Nne + F_Nne
    g = OmegaQ_Nne + F_Nne;
    % g % REMOVE
    if 1 == HBobj.isOsc
        % dOmegaQdX has an extra column for the derivative wrt to f, but dFdX does not; add a col of zeros
        dFdX(:, end+1) = sparse(size(dFdX,1),1);
    end
    J = dOmegaQdX + dFdX;
    evalsuccess = 1; % later, this should come up from device eval
end
% end of H_HB

function HBout = solve_HB(HBobj, Xinitguess_Nn, M, f0)
%function HBout = solve_HB(HBobj, Xinitguess_Nn, M, f0)
%- Xinitguess_Nn: the initial guess for NR, in Nn format
%  - if an oscillator, this should NOT include the frequency
%    guess as an unknown - that should be supplied as the f0
%    argument.
%- M: the number of 1-sided harmonics (excluding DC) to
%     be used in the 1-tone HB analysis.
%- f0: the fundamental frequency for non-autonomous HB analysis.
%      For autonomous systems (ie, isOsc == 1), f0 is used
%      as an initial guess for NR.
    HBobj.M = M;
    HBobj.N = 2*M+1;
    % set up U_Nni for HB analysis
    HBobj.U_Nni = setup_UNni(HBobj.DAE, f0, M, HBobj.N);

    HBobj.f0 = f0; 
    % fr is a guess if isOsc==1
    if 1 == HBobj.isOsc
        Xinitguess_Nn(end+1,1) = 1; % initial guess for fr, the relative frequency
    end

    % fixme: this printout should only happen at the appropriate dbglvl
    fprintf(2,'\n');
    [HBobj.solution, iters, HBobj.solvalid] = NR(@H_HB, [], Xinitguess_Nn, HBobj, HBobj.parms.NRparms);
    % TODO: have a HBparms and a HBparms.dbglvl
    if  HBobj.solvalid ~= 1
        fprintf(2,'\nHB solve: main NR failed after %d iterations\n', iters);
    else
        fprintf(2,'\nHB solve: main NR succeeded in %d iterations\n', iters);
    end

    %fprintf(2, 'JR, 2015/07/19: 2nd HB run to clean up complex conjugacy errors commented out\n');
    if 1==1 && 1 == HBobj.solvalid
        % fix complex conjugacy/realness errors that may have accumulated during complex NR
        % and run NR again. Hopefully we will get a solution that respects 
        % complex conjugacy/realness much more.

        n = feval(HBobj.DAE.nunks, HBobj.DAE);
        N = HBobj.N;
        if 1 == HBobj.isOsc
            sol_twoD = Nn_to_twoD(HBobj.solution(1:(end-1),1), n, N);
            fr = HBobj.solution(end,1);
            %f = HBobj.f0*fr;
        else
            sol_twoD = Nn_to_twoD(HBobj.solution, n, N);
        end

        M = HBobj.M;
        Xinitguess_twoD(:, 1) = real(sol_twoD(:,1)); % DC component
        avg = 0.5*(sol_twoD(:,2:(M+1)) + conj(sol_twoD(:,N:-1:(N-M+1))));
        Xinitguess_twoD(:, 2:(M+1)) = avg; % positive harmonics
        Xinitguess_twoD(:, N:-1:(N-M+1)) = conj(avg); % negative harmonics

        Xinitguess_Nn = twoD_to_Nn(Xinitguess_twoD);

        if 1 == HBobj.isOsc
            %Xinitguess_Nn(end+1,1) = real(fr);
            Xinitguess_Nn(end+1,1) = abs(fr);
            %Xinitguess_Nn(end+1,1) = real(fr)-imag(fr);
        end

        fprintf(2,'\nHB solve: re-running NR to try to clean up complex conjugacy/realness errors...\n', iters);
        [HBobj.solution, iters, HBobj.solvalid] = NR(@H_HB, [], Xinitguess_Nn, HBobj, HBobj.parms.NRparms);
        if  HBobj.solvalid ~= 1
            fprintf(2,'\nHB solve: clean-up NR failed after %d iterations\n', iters);
        else
            fprintf(2,'\nHB solve: clean-up NR succeeded in %d iterations\n', iters);
        end
    end

    [g, J, evalsuccess] = H_HB(HBobj.solution, HBobj);
    HBobj.residual_at_sol = g;
    HBobj.Jacobian_at_sol = J;
    if isfield(HBobj, 'S_HB')
        rmfield(HBobj, 'S_HB'); % previous sens. analysis no longer valid.
    end
    HBout = HBobj;
end % of solve_HB

function HBobj = compute_sensitivities_direct(HBobj, pObj)
%function HBobj = compute_sensitivities_direct(HBobj, pObj)
%- pObj - a Parameters object specifying parameters for the sensitivity 
%         analysis. See help Parameters.
%
%- computes the sensitivity matrix of the HB solution wrt the parameters and
%  stores it in HBobj.S_HB.
%  - see the derivation of the formula in
%    2017-03-19-Note-11-50--HB-parametric-adjoints.xoj
%
%- TODO - no thought has been paid to what the isOsc=1 case code really means
%         when writing this code, nor has it been tested for this case. It
%         probably won't work, or make any sense, for isOsc=1.
%
%
    % solution is stored in HBobj.solution by solve, but contains f if isOsc=1
    if 1 == HBobj.isOsc
        fr = HBobj.solution(end,1); % last entry of the unknown vector is the
                                    % relative frequency fr
        f = HBobj.f0*fr; 
        X_Nn = HBobj.solution(1:(end-1),1);
    else
        X_Nn = HBobj.solution;
        f = HBobj.f0; fr = 1;
    end
    % re-evaluate U (though this should have been set up before in HBobj.U_Nni)
    % TODO: there might be a bug in HBobj.solve for the isOsc=1 case. U needs
    % to change during the NR iteration since f changes, but does it? Check
    % this.
    U_Nni = setup_UNni(HBobj.DAE, f, HBobj.M, HBobj.N);

    J_HB = HBobj.Jacobian_at_sol; % stored by solve
    [dFdp, dQdp] = JpFQ_HB(X_Nn, U_Nni, pObj, HBobj); % derivatives of F and Q
                                                      % wrt parameters. Not
                                                      % isOsc=1 aware, but the
                                                      % F does not depend on f
                                                      % or fr, so this is fine.
    if 1 == HBobj.isOsc
        % dOmegaQdp has an extra column for the derivative wrt to f, but dFdX
        % does not; so, add a col of zeros to it
        dFdp(:, end+1) = sparse(size(dFdp,1),1);
    end

    [~, dOmegaQdp] = OmegaQ_HB([], dQdp, fr, HBobj); % re-use OmegaQ_HB to
                                                     % compute j 2 pi Omega_Nn
                                                     % dQdp. Does add a extra
                                                     % row at the end if
                                                     % isOsc=1.
    Jp_HB = dOmegaQdp + dFdp;

    HBobj.S_HB = - J_HB \ Jp_HB; % done: this is the direct sensitivity formula
end % of compute_sensitivities_direct

function hbsolStruct = getsolution_HB(HBobj)
%function hbsolStruct = getsolution_HB(HBobj)
%returns a structure containing the following fields:
% frequency: fundamental frequency of the HB analysis. If isOsc==1, this was solved for as part of HB.
% X_twoD: 2-D (n x N matrix) form of the HB solution - the Fourier coefficients in standard fft order.
% x_twoD: 2-D (n x N matrix) form of the HB solution, converted to the time domain (equispaced samples over the period).
% tpts: the timepoints corresponding to the columns of x_twoD
% U_twoD: 2-D (n x N matrix) form of the inputs to the DAE system - the Fourier coefficients in standard fft order.
% u_twoD: 2-D (n x N matrix) form of the inputs to the DAE system - time-domain samples.
% Jacobian_Nne_Nn: the Jacobian of the HB system H(X) (corresponding to Nn format for both H and X).
%   of size N*neqns x N*nunks if isOsc == 0,
%       N*neqns x (N*nunks+1) if isOsc == 1.
    if 1 ~= HBobj.solvalid
        fprintf(2,'HB getsolution: please run a successful solve first.\n')
        hbsolStruct = [];
        return;
    end
    n = feval(HBobj.DAE.nunks, HBobj.DAE);
    N = HBobj.N;
    if 0 == HBobj.isOsc
        hbsolStruct.frequency = HBobj.f0;
        hbsolStruct.f0 = HBobj.f0;
        hbsolStruct.X_twoD = Nn_to_twoD(HBobj.solution, n, N);
    elseif 1 == HBobj.isOsc
        fr = HBobj.solution(end,1);
        if is_not_small(imag(fr),real(fr), 1e-7, 1e-12)
            fprintf(2,'Warning: osc relative frequency fr has a significant imaginary component (%g/%g), ignoring\n', ...
                imag(fr), real(fr));
        end
        hbsolStruct.frequency = HBobj.f0*real(fr);
        hbsolStruct.fr = fr;
        hbsolStruct.f0 = HBobj.f0;
        hbsolStruct.X_twoD = Nn_to_twoD(HBobj.solution(1:(end-1),1), n, N);
    end
    hbsolStruct.residual_Nne = HBobj.residual_at_sol;
    hbsolStruct.Jacobian_Nne_Nn = HBobj.Jacobian_at_sol;
    hbsolStruct.tpts = 1/hbsolStruct.frequency*(0:(N-1))/N;
    xdash = N*ifft(full(hbsolStruct.X_twoD.')); % ifft works on cols
    hbsolStruct.x_twoD = real(xdash.');
    ni = feval(HBobj.DAE.ninputs, HBobj.DAE);
    if ni > 0
        hbsolStruct.U_twoD = Nn_to_twoD(HBobj.U_Nni, ni, N);
        udash = N*ifft(full(hbsolStruct.U_twoD.')); % ifft works on cols
        hbsolStruct.u_twoD = real(udash.');
    else
        hbsolStruct.U_twoD = sparse(0,N);
        hbsolStruct.u_twoD = sparse(0,N);
    end
    hbsolStruct.M = HBobj.M;
    hbsolStruct.N = N;
    hbsolStruct.n = n;
    hbsolStruct.ne = feval(HBobj.DAE.neqns, HBobj.DAE);
    hbsolStruct.ni = feval(HBobj.DAE.ninputs, HBobj.DAE);
end % of getsolution_HB

function S_HB = getSensitivities_HB(HBobj)
    if ~isfield(HBobj, 'S_HB')
        fprintf(1, 'HB.getSensitivities: you must run compute_sensitivities_direct() first.\n');
        S_HB = [];
        return;
    end
    S_HB = HBobj.S_HB;
end % getSensitivities_HB

function [figh, legends, colindex] = plot_HB(HBobj, stateoutputs, lgndprefix, linetype, figh, legends, colindex)
% Changelog: April 07, 2017
% Same change that was made to transientPlot recently. Until now, stateoutputs looked at DAE variables instead of DAE
% outputs for plotting. This has been changed. Instead of constructing a new C matrix from scratch, stateoutputs option
% now lets you pick from the outputs defined in the DAE, i.e., selecting a few rows from DAE.outputmatrix
	if 1 ~= HBobj.solvalid
		fprintf(2, sprintf('HB plot: please run solve successfully first.\n'));
		return;
	end

	DAE = HBobj.DAE;
	ninps = feval(DAE.ninputs, DAE);
	nunks = feval(DAE.nunks, DAE);
	if nargin < 2 || 0 == sum(size(stateoutputs))
		% plot DAE outputs
		C = feval(DAE.C, DAE);
		D = feval(DAE.D, DAE);
		onames = feval(DAE.outputnames, DAE);
	else % plot state outputs specified in stateoutputs
		% set up C, D, onames

                dae_C = DAE.C(DAE);
                dae_D = DAE.D(DAE);
                ninps = feval(DAE.ninputs, DAE);
                nunks = feval(DAE.nunks, DAE);

                varidxs = feval(stateoutputs.OutputIndices, stateoutputs);
                C = sparse(length(varidxs), nunks);
                for i=1:length(varidxs)
                        C(i,:) = dae_C(varidxs(i),:);
                end
		D = zeros(length(varidxs), ninps);
		onames = feval(stateoutputs.OutputNames, stateoutputs);
	end

	if nargin < 3 || 0 == sum(size(lgndprefix))
		lgndprefix = '';
	end

	if nargin < 4 || 0 == sum(size(linetype));
		linetype = '.-';
	end

	if nargin < 5 || 0 == sum(size(figh))
		figh = [];
	end

	if nargin < 6 || 0 == sum(size(legends))
		legends = {};;
	end

	if nargin < 7 || 0 == sum(size(colindex))
		colindex = 0;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (0 == strcmp('', lgndprefix))
		for i=1:length(onames)
			onames{i} = sprintf('%s: %s', lgndprefix, onames{i});
		end
	end

	%if length(legends) > 0
	%	onames = {legends{:}, onames{:}};
	%end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	sol = getsolution_HB(HBobj);
	noutputs = size(C,1); % ie, no of rows of C

	if 0 == noutputs
		fprintf(2,'HB plot: no outputs defined, nothing to do - returning.\n');
		return;
	end
	daename = feval(DAE.daename, DAE);
	M = HBobj.M;
	N = HBobj.N;

	out_FD_twoD = C*sol.X_twoD;
	if ninps > 0
		out_FD_twoD = out_FD_twoD + D*sol.U_twoD;
	end
	out_TD_twoD = C*sol.x_twoD;
	if ninps > 0
		out_TD_twoD = out_TD_twoD + D*sol.u_twoD;
	end
	out_TD_twoD(:,end+1) = out_TD_twoD(:,1); % repeat first point at end
						 % to make it visually periodic 
	tpts = 1/sol.frequency*(0:N)/N; % one more point than sol.tpts

	fharms = [0, 1:M, -M:-1];

	titlestr = sprintf('HB on %s:', feval(DAE.daename, DAE));
	time_units = DAE.time_units;
	[figh, legends, colindex] = plot_FD_TD(titlestr, onames, sol.frequency, fharms, out_FD_twoD, tpts, out_TD_twoD, time_units, figh, ...
            legends, colindex);
end
% end of plot_HB

function U_Nni = setup_UNni(DAE, f, M, N)
%function u_Nni = setup_UNni(DAE, f0, M, N)
%set up U_Nn for HB analysis, checking the sizes
%of the DAE's HB inputs against N and adjusting
%as necessary.
	% set up U_Nni for HB
	ninputs = feval(DAE.ninputs, DAE);
	if ninputs > 0
		Us = feval(DAE.uHB, f, M, DAE);
		if size(Us, 1) ~= ninputs
			error(sprintf('DAE''s HB input vector size %d != no of DAE inputs %d', size(Us,1), ninputs));
		end
		inN = size(Us,2);

		% check that the number of harmonics in the DAE's input vector is odd
		if 1 ~= mod(inN,2) % not an odd number of harmonics
			error(sprintf('DAE''s HB input vector no of total harmonics (%d) not odd', inN));
		end

		inM = (inN-1)/2;

		% adjust the number of harmonics of Us to be exactly N
		if inN < N
			nonconjs = Us(:,1:(inM+1));
			conjs = Us(:,(inM+2):end);
			fills = sparse(ninputs,2*(M-inM));
			Us = [nonconjs, fills, conjs]; % no of cols should be exactly N, now
		elseif inN > N
			fprintf(2, 'Warning: HB using fewer harmonics (%d) than DAE''s HB input (%d) - extra input harmonics ignored',...
				HBobj.M, inM);
			nonconjs = Us(:,1:(M+1));
			conjs = Us(:,(end-M+1):end)
			Us = [nonconjs, conjs]; % no of cols should be exactly N, now
		end

		% store the HB inputs
		U_Nni = twoD_to_Nn(Us);
	else
		U_Nni = [];
	end
end
% end of setup_UNn
