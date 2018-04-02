function [dFdp, dQdp] = JpFQ_HB(X_Nn, U_Nni, pObj, HBobj)
%function [dFdp, dQdp] = JpFQ_HB(X_Nn, U_Nni, pObj, HBobj)
%
%This function provides parametic sensitivity matrices for the HB equation
%H_HB(X, U, p) = 0. More precisely, we have
%   H_HB(X, U, p) = j*2*pi*f0*Omega*Q(X, p) + F(X, p);
%this function returns dFdP = d/dp F(X,p), and dQdp = d/dp Q(X,p).
%For details, see 2017-03-19-Note-11-50--HB-parametric-adjoints.xoj.
%
%The arguments X and U should be in Nn format (unlike in the derivation, which
%uses nN format). pObj should be a Parameters object (see help Parameters)
%containing the parameters with respect to which derivatives are desired.
%
%The outputs dFdp and dQdp are matrices matrix of size Nn rows x ni columns.
%The Nn rows are in Nn format; the columns are in the same order as the
%parameters in pObj.
%
%recall:
%	nN format: N block-vectors (of F coeffs in standard order, or uniformly
%	           space time samples) stacked one after the other. Each
%	           block-vector is of size n (in the order of the DAE's unknowns)
%	Nn format: n block-vectors (in the order of the DAE's unknowns) stacked one
%	           after the other.  Each block-vector is of size N (the F coeffs
%	           in standard order, or uniformly spaced time samples) of the
%	           corresponding DAE unknown.
%
%For a derivation of the equations implemented here, see
%2017-03-19-Note-11-50--HB-parametric-adjoints.xoj.
%
%TODO: most of the code in this is exactly repeated from FQJ_HB.m. Should
%re-organize both files to prevent such repetition.
%
%Examples
%--------
%
%TODO
%
%
%See also
%--------
%
%TODO

%Author: Jaijeet Roychowdhury <jr@berkeley.edu>, 2017/03/19. Copied FQJ_HB.m
%        and modified.
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2008-2017 Jaijeet Roychowdhury <jr@berkeley.edu>. Some rights
%               reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	DAE = HBobj.DAE;
	N = HBobj.N; % no of harmonics
	n = feval(DAE.nunks, DAE); % no of unknowns
	ni = feval(DAE.ninputs, DAE); % no of inputs
	ne = feval(DAE.neqns, DAE); % no of equations (size of DAE.f(), DAE.q())


	%%%%%%%%%% step 1: X_Nn to x_Nn, and U_Nn to u_Nn, using the ifft

	X_twoD = full(Nn_to_twoD(X_Nn, n, N));
	x_twoDdash = N*ifft(X_twoD.'); % ifft works on each column
	x_Nn = twoD_to_Nn(x_twoDdash.');

	U_twoD = full(Nn_to_twoD(U_Nni, ni, N));
	u_twoDdash = N*ifft(U_twoD.'); % ifft works on each column
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

	%%%%%%%%%% step 3: compute time-domain DAE parameter Jacs at each timepoint
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

		% Jacobians: dfdp, dqdp at each timepoint
		if 1 == DAE.f_takes_inputs
            [dfdp_ti, dqdp_ti] = feval(DAE.dfq_dp, x, u, pObj, 'fq', DAE);
		else
            [dfdp_ti, dqdp_ti] = feval(DAE.dfq_dp, x, pObj, 'fq', DAE);
		end

        % stack up the Jacobians
		dfdp(startidx_neN:endidx_neN,:) = dfdp_ti;
		dqdp(startidx_neN:endidx_neN,:) = dqdp_ti;
	end
	%%%%%%%%%% end step 3

	%%%%%%%%%% step 4: permute rows from neN to Nne format

	% set up permutation matrix to_Nne_P for the Jacobians 
    % hit with to_Nne_P from the left to permute the rows (later)
	to_Nne_indices = to_Nn(1:ne*N,ne,N);
	eye_neN = speye(ne*N);
	to_Nne_P = eye_neN(to_Nne_indices,:);
	%%%%%%%%%% end step 4

	%%%%%%%%%% step 5: f_Nne to F_Nne using the fft
	% Jacobian: dFdf_Nne, dQdq_Nne are just n blocks of NxN DFT matrices
	DFT_N = fft(eye(N));
	dFQdfq_Nne = sparse(ne*N,ne*N);
	for i=1:ne
		startidx_Nne = (i-1)*N+1;
		endidx_Nne = i*N;
		dFQdfq_Nne(startidx_Nne:endidx_Nne,startidx_Nne:endidx_Nne) = 1/N*DFT_N;
	end
	%%%%%%%%%% end step 5

	% build the end-to-end Jacobians
	dFdp = dFQdfq_Nne * to_Nne_P * dfdp;
	dQdp = dFQdfq_Nne * to_Nne_P * dqdp;
	% for iterative methods, each of these matrices can be applied quickly to a vector
end
