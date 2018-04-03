function [x_dct, T] = m_DCT2(x,mode)
% function [x_dct, T] = m_DCT2(x,mode)
% Author: Archit Gupta (September 06, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Implementation of DCT-II. DCT transform of the second kind.
%   Refer WIKIPEDIA or <CS_DIR>/4_Notes/Discrete-cosine-transform.xoj for more
%   details. For a given sequence of values x[n], this computes the cosine
%   transform of an even extension about points halfway from the left and right
%   ends.
%   
%   INPUTS:
%       x - A sequence of values for which DCT-2 has to be computed [Expect a
%       column vector. A row vector will produce an ERROR]
%       mode (optional) - This can be left empty if the regular DCT is to be
%       computed. However, if we want to compute the DCT matrix by using the DFT
%       Available options: 'MODE_EXTENDED', 'CHEB_SERIES'
%
%   OUTPUTS:
%       x_dct - A sequence of DCT coefficients (same length as x)
%       T (optional) - The tranform matrix associated with the size of input x
%
%   NOTE: This code is experimental. It does not attain the complexity of
%   O{nlog(n)}, where n is the length of the vector x. This has a complecity of
%   O(n^2), which is probably good enough for the small scale experiments that
%   we are doing. For better efficiency (atleast better asymptotic efficiency)
%   use DFT instead with an appropriately extended sequence.
%
%   TRANSFORM:
%       X(k) = sum_{n=0}^{N-1} x(n)*cos(pi*(n+0.5)*k)/N
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    s_x = size(x);
    if ( length(s_x) > 2 )
        error(['DCT for inputs having more than 2 dimensions', ...
            ' not supported at the moment']);
        x_dct = 0;
        T = 0;
        return;
    end


    if (s_x(1) == 1)
        fprintf(2, ['WARNING: If row vector was given the output is wrong.\n', ...
            'Otherwise DCT2 is trivial, what is wrong with you!\n']);
        x_dct = x;
        T = 1;
        return;
    end

    % Baseline code: MATLAB implementation (and several others) require scaling of
    % the overall matrix and also the X0 term (so that the transform matrix is
    % orthogonal). See scaled implemetation here for details. The version below
    % directly corresponds to a DFT of a sequence of length 4N where an even
    % extension at a halfway point is placed at the even point in the sequence.

    n_elems = s_x(1);
    % Let k denote the index of the DCT coefficient and let n denote the index of
    % the time domain coefficient
    k = [0:n_elems-1]';
    n = [0:n_elems-1]';

    T = cos(pi*k*(n+0.5)'/n_elems);

    % Standard - orthogonal DCT transform
    % Scaled implementation - Usual DCT
    T = sqrt(2/n_elems)*T;
    T(1,:) = T(1,:)/sqrt(2);

    x_dct = T * x;
end