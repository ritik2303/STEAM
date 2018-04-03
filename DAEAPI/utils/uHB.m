function uout = uHB(f, M, DAE)
% function uout = uHB(f, DAE)
%Author: Jaijeet Roychowdhury <jr@berkeley.edu> 2012/06/08
%Changelog: April 06, 2017
%   - Currently, calls to this function do not include the number of harmonics that are expected. So, it is not possible
%   to fill in zeros in the proper FFT order. This is particularly problematic if there are multiple inputs and only 1
%   of them is a non-DC HB input   
%   > This has been changed. Calls to uHB now need to supply M, the number of Harmonics so that the uout matrix can be
%   appropriately sized when there are multiple inputs. For example: A for circuit with VDD and VIN, VDD has only 1
%   harmonic (non-zero), specified by the function associated with its set_uHB byt VIN has atleast 3 harmonics.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Type "help MAPPlicense" at the MATLAB/Octave prompt to see the license      %
%% for this software.                                                          %
%% Copyright (C) 2008-2013 Jaijeet Roychowdhury <jr@berkeley.edu>. All rights  %
%% reserved.                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	nI = length(DAE.uHBfunc_updates.indices);

	if isempty(DAE.uHBfunc_updates.vector_uHBfunc)
		uout = zeros(DAE.ninputs(DAE), 2*M + 1);
	else
		uout = feval(DAE.uHBfunc_updates.vector_uHBfunc, f, DAE.uHBfunc_updates.vector_uHBfunc_args);
	end

	for i=1:nI
		u = feval(DAE.uHBfunc_updates.uHBfunclist{i}, f, DAE.uHBfunc_updates.uHBargslist{i});
                n_coeffs_obtained = length(u);
                u_full = zeros(1, 2*M+1);
                if ( rem(n_coeffs_obtained,2) == 0 )
                    error('uHB: HB Function returned a even number of Fourier Coefficients');
                end
                n_non_dc_coeffs = floor(n_coeffs_obtained/2); % Using a floor to ensure there are no rounding issues
                % The other option is to use (x-1)/2 but it is not an integer categorically

                u_full(1,1) = u(1);
                u_full(1,2:1+n_non_dc_coeffs) = u(2:1+n_non_dc_coeffs);
                u_full(1,end-n_non_dc_coeffs+1:end) = u(n_non_dc_coeffs+2:end);
		% u should be a row vector of F coeffs
		uout(DAE.uHBfunc_updates.indices(i),:) = u_full;
		% FIXME: we want to resize each u to the max row length, filling in zeros in the proper FFT order
	end
end % uHB