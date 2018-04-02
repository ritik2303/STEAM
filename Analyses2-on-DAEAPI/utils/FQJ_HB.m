function [F_Nne, Q_Nne, dFdX, dQdX] = FQJ_HB(X_Nn, U_Nni, HBobj)
%function [F_Nne, Q_Nne, dFdX, dQdX] = FQJ_HB(X_Nn, U_Nni, HBobj)
%This takes in X and U in Nn format and produces F/Q/dFdX/dQdX in the same format.
%recall:
%	nN format: N block-vectors (of F coeffs in standard order, or uniformly space time 
%	      samples) stacked one after the other. Each block-vector is of size 
%	      n (in the order of the DAE's unknowns)
%	Nn format: n block-vectors (in the order of the DAE's unknowns) stacked one after the other. 
%	      Each block-vector is of size N (the F coeffs in standard order, or uniformly
%	      spaced time samples) of the corresponding DAE unknown.
%The Jacobians are block-sparse - the sparsity pattern of DAE's Jacobians is retained at the
%block level; each scalar point becomes a dense NxN circulant matrix
%
%Author: Jaijeet Roychowdhury <jr@berkeley.edu> 2012/06/04
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J. Roychowdhury.
% Copyright (C) 2008-2012 Jaijeet Roychowdhury <jr@berkeley.edu>. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% global HBdebuglvl; % REMOVE

	DAE = HBobj.DAE;
	N = HBobj.N; % no of harmonics
	n = feval(DAE.nunks, DAE); % no of unknowns
	ni = feval(DAE.ninputs, DAE); % no of inputs
	ne = feval(DAE.neqns, DAE); % no of equations (size of DAE.f(), DAE.q())


	%%%%%%%%%% step 1: X_Nn to x_Nn, and U_Nn to u_Nn, using the ifft

	X_twoD = full(Nn_to_twoD(X_Nn, n, N));
	x_twoDdash = N*ifft(full(X_twoD).'); % ifft works on each column
	x_Nn = twoD_to_Nn(x_twoDdash.');

	%{
	if 1 == HBdebuglvl 
		x_twoD = x_twoDdash.'
	end
	%}


	U_twoD = full(Nn_to_twoD(U_Nni, ni, N));
	u_twoDdash = N*ifft(full(U_twoD).'); % ifft works on each column
	u_Nni = twoD_to_Nn(u_twoDdash.');

	% alternative form below does the above using a for loop, 
	%	useful for understanding, and writing Jacobian code
	%% first, convert X_Nn to x_Nn using the ifft
	%for j=1:n
	%	startidx_Nn = (j-1)*N+1;
	%	endidx_Nn = j*N;
	%	x_Nn(startidx_Nn:endidx_Nn) = N*ifft(X_Nn(startidx_Nn:endidx_Nn));
	%	% above could be done in place, overwriting X_Nn, if memory is an issue
	%end
	%% now, x_Nn contains the time-domain samples of x(t) in Nn format

	% dxdX_Nn is just n blocks of NxN IDFT matrices
	IDFT_N = ifft(eye(N));
	dxdX_Nn = sparse(n*N,n*N);
	for i=1:n
		startidx_Nn = (i-1)*N+1;
		endidx_Nn = i*N;
		dxdX_Nn(startidx_Nn:endidx_Nn,startidx_Nn:endidx_Nn) = N*IDFT_N;
	end

	%%%%%%%%%% end step 1

	%%%%%%%%%% step 2: permute x_Nn to x_nN, u_Nni to u_niN
	x_nN = to_nN(x_Nn, n, N); % equivalent to: to_nN_P * x_Nn;
	u_niN = to_nN(u_Nni, ni, N); % equivalent to: to_niN_P * u_Nni;

	% Jacobian: dxnN_dxNn = to_nN_P
	to_nN_indices = to_nN(1:n*N,n,N);
	eye_nN = speye(n*N);
	to_nN_P = eye_nN(to_nN_indices,:);
	%%%%%%%%%% end step 2

	%%%%%%%%%% step 3: compute f/q_neN by evaluating DAE.f/q() at each timepoint
        %{
        flag.f =1; flag.q =1; 
        flag.dfdx =1; flag.dqdx =1;
        % flag.dfdu = 1; % JR 2016/03/01 - took this out, not used in LMS!
        flag.dfdu = 0;
        if 1 == DAE.f_takes_inputs
            % flag.dfdu = 1; % JR 2016/03/01 - not used in LMS!
            fqJout = feval(DAE.fqJ, x, unew, flag, DAE);
            fout = fqJout.f;
            % dfdu = fqJout.dfdu; % JR 2016/03/01 - not used in LMS!
        else
            fqJout = feval(DAE.fqJ, x, flag, DAE);
            Bmat = feval(DAE.B, DAE);
            if ~isempty(Bmat)
                fout = fqJout.f + Bmat * unew;
            else
                fout = fqJout.f;
            end
            % dfdu = Bmat; % JR 2016/03/01 - not used in LMS!
        end
        qout = fqJout.q;
        dfdx = fqJout.dfdx;
        dqdx = fqJout.dqdx;
        %}

    % set up flags for DAE.fqJ call
    flag.f =1; flag.q =1; 
    flag.dfdx =1; flag.dqdx =1;
    flag.dfdu = 0; % not needed for HB NR solution

	if 0 == DAE.f_takes_inputs
		DAE_B = feval(DAE.B, DAE);
	end

	for i=1:N
		startidx_nN = (i-1)*n+1; endidx_nN = i*n;
		startidx_niN = (i-1)*ni+1; endidx_niN = i*ni;
		x = x_nN(startidx_nN:endidx_nN,1);
		u = u_niN(startidx_niN:endidx_niN,1);

		% make x and u real; check that their imaginary components are negligible
		reltol = 1e-9; abstol = 1e-12;
		if is_not_small(imag(x),real(x),reltol,abstol)
			fprintf(2,'HB warning: x(t_%d) has significant imaginary components (%g/%g) - ignoring them.\n', ...
				i, norm(imag(x)), norm(real(x)));
		end
		x = real(x);

		if length(u) > 0 
			if is_not_small(imag(u),real(u),reltol,abstol)
				fprintf(2,'HB warning: u(t_%d) has significant imaginary components (%g/%g) - ignoring them.\n', ...
				i, norm(imag(u)), norm(real(u)));
			end
			u = real(u);
		end

		startidx_neN = (i-1)*ne+1; endidx_neN = i*ne;

        %{
        OLD CODE, separate calls to f/q/Jf/Jq, inefficient; commented out 2017/04/07, JR.
		if 1 == DAE.f_takes_inputs
			f_neN(startidx_neN:endidx_neN,1) = feval(DAE.f, x, u, DAE);
			
			%if 1 == HBdebuglvl
			%	xi = x
			%	ui = u
			%	f_xi = feval(DAE.f, x, u, DAE)
			%end
			
		else
			f_neN(startidx_neN:endidx_neN,1) = feval(DAE.f, x, DAE) + DAE_B*u;
		end
		q_neN(startidx_neN:endidx_neN,1) = feval(DAE.q, x, DAE);

		% Jacobians: dfdx, dqdx
		if 1 == DAE.f_takes_inputs
			dfdx(startidx_neN:endidx_neN,startidx_nN:endidx_nN) = feval(DAE.df_dx, x, u, DAE);
		else
			dfdx(startidx_neN:endidx_neN,startidx_nN:endidx_nN) = feval(DAE.df_dx, x, DAE);
		end
		dqdx(startidx_neN:endidx_neN,startidx_nN:endidx_nN) = feval(DAE.dq_dx, x, DAE);
        END OLD CODE
        %}
        if 1 == DAE.f_takes_inputs
            fqJout = feval(DAE.fqJ, x, u, flag, DAE);
            fout = fqJout.f;
        else
            fqJout = feval(DAE.fqJ, x, flag, DAE);
            if ~isempty(Bmat)
                fout = fqJout.f + DAE_B * u;
            else
                fout = fqJout.f;
            end
        end

		f_neN(startidx_neN:endidx_neN,1) = fout;
		q_neN(startidx_neN:endidx_neN,1) = fqJout.q;
		dfdx(startidx_neN:endidx_neN,startidx_nN:endidx_nN) = fqJout.dfdx;
		dqdx(startidx_neN:endidx_neN,startidx_nN:endidx_nN) = fqJout.dqdx;
	end
	%%%%%%%%%% end step 3

	%%%%%%%%%% step 4: permute f_neN to f_Nne
	f_Nne = to_Nn(f_neN, ne, N); % equivalent to: to_Nne_P * f_neN;
	q_Nne = to_Nn(q_neN, ne, N); % equivalent to: to_Nne_P * q_neN;

	%{
	if 1 == HBdebuglvl
		f_twoD = nN_to_twoD(f_neN, ne, N)
		q_twoD = nN_to_twoD(q_neN, ne, N)
	end
	%}

	% Jacobian: df/qNne_dfneN = to_Nne_P
	to_Nne_indices = to_Nn(1:ne*N,ne,N);
	eye_neN = speye(ne*N);
	to_Nne_P = eye_neN(to_Nne_indices,:);
	%%%%%%%%%% end step 4

	%%%%%%%%%% step 5: f_Nne to F_Nne using the fft
	f_twoD = Nn_to_twoD(f_Nne, ne, N);
	% f_twoD % REMOVE
	F_twoDdash = 1/N*fft(full(f_twoD).'); % fft works on each column
	F_Nne = twoD_to_Nn(F_twoDdash.');

	q_twoD = Nn_to_twoD(q_Nne, ne, N);
	Q_twoDdash = 1/N*fft(full(q_twoD).'); % fft works on each column
	Q_Nne = twoD_to_Nn(Q_twoDdash.');

	% Jacobian: dFdf_Nne, dQdf_Nne are just n blocks of NxN DFT matrices
	DFT_N = fft(eye(N));
	dFQdfq_Nne = sparse(ne*N,ne*N);
	for i=1:ne
		startidx_Nne = (i-1)*N+1;
		endidx_Nne = i*N;
		dFQdfq_Nne(startidx_Nne:endidx_Nne,startidx_Nne:endidx_Nne) = 1/N*DFT_N;
	end
	%%%%%%%%%% end step 5

	% build the end-to-end Jacobians
	dFdX = dFQdfq_Nne * to_Nne_P * dfdx * to_nN_P * dxdX_Nn;
	dQdX = dFQdfq_Nne * to_Nne_P * dqdx * to_nN_P * dxdX_Nn;
	% for iterative methods, each of these matrices can be applied quickly to a vector
end
