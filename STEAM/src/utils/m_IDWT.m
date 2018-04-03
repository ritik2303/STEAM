function x_out = m_IDWT(ca, cd, basis)
% function x_out = m_IDWT(ca, cd, basis) or
% function x_out = m_IDWT(cacd, basis)
% Author: Archit Gupta (September 26, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementation of the inverse discrete wavelet transform. Currently, we are
% only dealing with a single level Haar wavelet, which is indeed a tight frame.
% So, the reconstruction is rather trivial. If the approximation coefficients
% and the detail coefficients are provided, we can easily reconstruct the
% original vector x using the formula:
%                       x = A^(-1) <x,v_{i}>v_{i}
% For the case of the single-level Haar wavelets, the frame bound is sqrt(1/2).
% INPUTS:
%   ca: Approximation coefficients
%   cd: Detail coefficients
%   basis: Currently 'haar' is the only accepted basis
%
% ALTERNATE USAGE: Keeping in mind that this would be used along with compressed
% sensing, the user can choose to provide the approximation and the detail
% coefficients in a single vector. The EXPECTED ORDER IS [CA; CD]. Both ca and
% cd should be column vectors in the first place. 
%
% If only one input is provided, the "default" basis will be assumed and it will
% also be assumed that a concantenation of the approximation and detail
% coefficients has been provided.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Setting up defaults
    def_basis = 'haar';
    
    nCoeffs = length(ca);
    % Processing input arguments
    if (nargin < 3)
        if (nargin < 2)
            basis = def_basis;
        else
            basis = cd;
        end
    
        nCoeffs = length(ca)/2;
        % If you are seeing an error from one of the two following lines of code,
        % double-check the size of the concatenated vector that is being input. It
        % might be an odd number, or worse, 0!
        cd = ca(1+nCoeffs:end,1);
        ca = ca(1:nCoeffs,1);
    end

    if (~strcmp(basis, 'haar'))
        fprintf(2, ['Only the Haar Basis is supported at the moment\n.' ...
            'See help m_DWT for details\n']);
        x_out = [];
        return;
    end

    M = 2*nCoeffs;
    x_out = zeros(M,1);
    index = [0:M-1]';
    j0 = 0;

    for k = 0:2:M-1
        x_out = x_out + ca(1+k/2,1)*haar_phi(index, j0, k)/sqrt(2);
        x_out = x_out + cd(1+k/2,1)*haar_psi(index, j0, k)/sqrt(2);
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
