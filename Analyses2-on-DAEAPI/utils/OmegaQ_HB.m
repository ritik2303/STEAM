function [OmegaQ_Nne, dOmegaQdX] = OmegaQ_HB(Q_Nne, dQdX, fr, HBobj)
%function [OmegaQ_Nne, dOmegaQdX] = OmegaQ_HB(Q_Nne, dQdX, fr, HBobj)
%
%This function returns OmegaQ, ie, d/dt q(x(t)) in terms of Fourier coeffs in Nn
%vector format, given the Fourier coeffs of q(x(t)) in Nn format.  It also
%returns the Jacobian of OmegaQ with respect to X_Nn, (ie x(t) in Nn format)
%given dQdX.
%


%% Author: Jaijeet Roychowdhury <jr@berkeley.edu> 2012/06/05

%% Changelog:

% 2017/03/19: small update to support empty first arg for Jacobian only
%            processing. JR.


% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2008-2017 Jaijeet Roychowdhury <jr@berkeley.edu>. Some rights
% reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	ne = feval(HBobj.DAE.neqns, HBobj.DAE);

    if ~isempty(Q_Nne) % Q_Nne can be [] if only the Jacobian is needed
        Q_twoD = Nn_to_twoD(Q_Nne, ne, HBobj.N);
        OmegaCoeffs = (0+1i)*2*pi*[0,1:HBobj.M,-HBobj.M:-1]; % d/dt operator on Fourier coeffs; row vector
        OmegaQ_twoD = (ones(size(Q_twoD,1),1) * OmegaCoeffs) .* Q_twoD;
        OmegaQ_Nne = HBobj.f0*fr*twoD_to_Nn(OmegaQ_twoD);
    else
        OmegaQ_Nne = [];
    end

	%{
	% a more explicit way of doing the above, using a for loop
	Q_neN = to_nN(Q_Nne, ne, HBobj.N);
	for j=1:HBobj.N
		if j <= HBobj.M+1
			F_idx = j-1;
		else
			F_idx = -(HBobj.N-j+1);
		end
		startidx_neN = (j-1)*ne+1:
		endidx_neN = j*ne;
		OmegaQ_neN(startidx_neN:endidx_neN,1) = (0+1i)*2*pi*HBobj.f0*fr*F_idx*Q_neN(startidx_neN:endidx_neN,1);
	end
	OmegaQ_Nne = to_Nn(OmegaQ_neN, ne, HBobj.N);
	%}


	% Jacobian:
	%
	% from the above explicit way, it is clearl that: dOmegaQ_dX = to_Nne_P * dOmegaQ_dQ_neN * to_neN_P * dQ_dX
	% where dOmegaQ_dQ_neN = jay*2*pi*f*blkdiag( [0,0,..,0], [1,1,..,1], ..., [M,M,..,M], [-M,-M,..,-M], ..., [-1,-1,..,-1])
	% dOmegaQ_dQ_neN can be generated conveniently using kron, as follows:
	%	Fmults = (0+1i)*2*pi*HBobj.f0*fr*sparse(diag(sparse([0,1:M,-M:-1])))
	%	eye_ne = speye(ne);
	%	dOmegaQ_dQ_neN = kron(Fmults, eye_ne); % the order Fmults, eye_ne is important
	%
	% Indeed, we can directly get to_Nne_P * dOmegaQ_dQ_neN * to_neN_P simply by reversing the kron arg order:
	%	to_Nne_P * dOmegaQ_dQ_neN * to_neN_P = kron(eye_ne, Fmults)
	%
	Fmults = (0+1i)*2*pi*HBobj.f0*fr*sparse(diag(sparse([0,1:HBobj.M,-HBobj.M:-1])));
	eye_ne = speye(ne);
	dOmegaQdX = kron(eye_ne, Fmults) * dQdX;

	if 1 == HBobj.isOsc % oscillator
		% fr is an unknown - add derivatives wrt fr as last col of dOmegaQdX
		% (last col because by convention, the freq unknown is added to
		%  the end of X_Nn)
		% from the code for OmegaQ_Nne above, we have:
		%   	dOmegaQ_df = HBobj.f0*twoD_to_Nn(OmegaQ_twoD)
		dOmegaQdX(:,end+1) = HBobj.f0*twoD_to_Nn(OmegaQ_twoD);
	end
end
