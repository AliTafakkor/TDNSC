function Z = squareform(Y,dir) %#codegen
%SQUAREFORM Reformat a distance matrix between upper triangular and square form.
%   Z = SQUAREFORM(Y), if Y is a vector as created by the PDIST function,
%   converts Y into a symmetric, square format, so that Z(i,j) denotes the
%   distance between the i and j objects in the original data.
%
%   Y = SQUAREFORM(Z), if Z is a symmetric, square matrix with zeros along
%   the diagonal, creates a vector Y containing the Z elements below the
%   diagonal.  Y has the same format as the output from the PDIST function.
%
%   Z = SQUAREFORM(Y,'tovector') forces SQUAREFORM to treat Y as a matrix.
%   Y = SQUAREFORM(Z,'tomatrix') forces SQUAREFORM to treat Z as a vector.
%   These formats are useful if the input has a single element, so it is
%   ambiguous as to whether it is a vector or square matrix.
%
%   Example:  If Y = (1:6) and X = [0  1  2  3
%                                   1  0  4  5
%                                   2  4  0  6
%                                   3  5  6  0],
%             then squareform(Y) is X, and squareform(X) is Y.
%
%   See also PDIST.

%   Copyright 2017 The MathWorks, Inc.
    coder.internal.prefer_const
    coder.internal.errorIf(~(isnumeric(Y) || islogical(Y)) || ~ismatrix(Y),'stats:squareform:BadInput');

    if nargin<2 || isempty(dir)
        if isvector(Y)
            direction = MATRIX;
        else
            direction = VECTOR;
        end
    else
        coder.internal.errorIf(~coder.internal.isConst(dir),'stats:coder:pdist2:ExpectedConstant','Direction');
        validateattributes(dir,{'string','char'},{'row'},mfilename,'DIRECTION');
        dirtemp = validatestring(dir,{'tomatrix','tovector'},mfilename,'DIRECTION');
        if strncmpi(dirtemp,'tov',3)
            direction = VECTOR;
        else
            direction = MATRIX;
        end
    end

    [m, n] = size(Y);

    switch direction
      case MATRIX
        coder.internal.errorIf(~isvector(Y),'stats:squareform:BadInputVector');

        if m~=1
            n = m;
        end

        m = ceil(sqrt(2*n)); % this is an approximation to the actual solution of
                             % the quadratic equation (m choose 2) = n i.e. (1 + sqrt(1+8*n))/2.
        coder.internal.errorIf( m*(m-1)/2 ~= n,'stats:squareform:BadVectorSize');

        Z = zeros(m,'like',Y);

        if m>1
            k = 1;
            for j = 1:m-1
                for i = j+1:m
                    % The ctranspose on a scalar below is for conjugate
                    % in the complex case, and unlike CONJ, it is
                    % safe with logicals.
                    Z(i,j) = Y(k);
                    Z(j,i) = Y(k)';
                    k = k + 1;
                end
            end
        end

      otherwise % case VECTOR
        coder.internal.errorIf( m~=n || ~all(diag(Y)==0),'stats:squareform:BadInputMatrix');
        Z = zeros(1,n*(n-1)/2, 'like',Y);
        k = 1;

        for j = 1:m-1
            for i = j+1:m                % here conjugate operation is necessary since the lower
                Z(k) = Y(i,j)';  % triangle is being traversed (for speed) and the upper
                k = k + 1;       % triangle needs to be returned
            end
        end
    end
end

function out = VECTOR
    out = uint8(2);
end

function out = MATRIX
    out = uint8(1);
end
