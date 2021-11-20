function [ strideSigma ] = GetStride( k, j, n, p, m )

% [ strideSigma ] = GetStride( k, j, n, p, m )
% Purpose: Returns the stride to access Sigma^{(k)}_{j} in the 1D array Sigma.
% Created:     09.11.2016
% Last change: 09.11.2016

if or(k<=0,k>m)
    error('k must be 1<=k<=m.');
end

if or(j==1,j==2)
    a = j-1;
    strideSigma = (2*(k-1)+a)*n*p+1:(2*k-1+a)*n*p;
else
    error('j must be either 1 or 2.');
end

end