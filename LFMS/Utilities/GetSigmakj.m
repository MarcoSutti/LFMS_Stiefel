function [ reshapedSigmakj ] = GetSigmakj( Sigma, k, j, n, p, m )

% [ reshapedSigmakj ] = GetSigmakj( Sigma, k, j, n, p, m )
% Purpose: Returns Sigma^{(k)}_{j} in reshaped form.
%          k is the index of the subinterval.
%          j=1, Stiefel point; j=2, tangent vector.
% Created:     09.11.2016
% Last change: 09.11.2016

if or(k<=0,k>m)
    error('k must be 1<=k<=m.');
end

if or(j==1,j==2)
    reshapedSigmakj = reshape( Sigma( GetStride( k, j, n, p, m ) ), [n,p] );
else
    error('j must be either 1 or 2.');
end

end