function [ca, cd] = m_DWT(x_in, basis)
% function [ca, cd] = m_DWT(x_in, basis)
% Author: Archit Gupta (September 26, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementation of the discrete wavelet transform. Please refer to
% <CS_DIR>/4_Notes/Discrete_Wavelet_Transform.xoj for details of implementation.
% Or look at dissecting_haar_wavelet for details on the Haar wavelet
% INPUTS:
%   x_in: Input for which the discrete wavelet transform must be calculated
%   basis: The name of the basis which should be used for computing the wavelet
%   transform.  Currently, only the 'haar' basis is supported
%
% OUTPUTS:
%   ca: Approximation coefficients - Convolution of x_in with the "low-pass"
%   filter corresponding to the transform basis. Only the even coefficients are
%   returned
%   cd: Detail coefficients - Convolution of x_in with the "high-pass" filter
%   corresponding to the transform basis. Again, only the even coefficients are
%   returned.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (~strcmp(basis, 'haar'))
        fprintf(2, ['Only the Haar Basis is supported at the moment\n.' ...
            'See help m_DWT for details\n']);
        ca = [];
        cd = [];
        return;
    end

    M = length(x_in);
    index = [0:M-1]';
    nCoeffs = ceil(M/2);
    ca = zeros(nCoeffs,1);
    cd = zeros(nCoeffs,1);
    j0 = 0;         % We are computing a single-level wavelet transform, so j0
                    % takes just one value, i.e. 0
    for k = 0:2:M-1
        ca(1+k/2,1) = dot(x_in, haar_phi(index, j0, k))/sqrt(2); 
        cd(1+k/2,1) = dot(x_in, haar_psi(index, j0, k))/sqrt(2); 
    end

    if (nargout < 2)
        ca = [ca; cd];  % return concatenated coefficients. Useful for CS
    end
end

% Some basic function definitions required for computing the haar wavelets

function dotproduct = dot(vec_x, vec_y)
% Returns the sum of the product of pairwise elements
    dotproduct = sum(vec_x .* vec_y);
end

function y = haar_phi(x,j,k)
% phi(x) is a symmetric rectangle centered at 1/2 and has a width 1. The
% function phi_{j,k}(x) = 2^(j/2) * phi( ((2^j) * x) - k), which is a dilation
% by a factor of 2^j and a translation by k units. In discrete domain, where x
% can only take integer values, however, the function remains the same.
    x_new = (2^j)*x - k;
    y = (2^(j/2)) * (x_new >= 0) .* (x_new <= 1);
end

function y = haar_psi(x,j,k)
% psi(x) is a anti-symmetric function centered at 1/2. Its value is 1 for 0 <=
% x < 1/2 and -1 for 1/2 <= x <= 1. psi_{j,k}(x) is defined in terms of psi(x) in
% the  same way as phi(x).
    x_new = (2^j)*x - k;
    y = (2^(j/2)) * (x_new >= 0).* ( (x_new < 0.5) + (-1)*(x_new>=0.5) ) .* ...
        (x_new <= 1);
end
