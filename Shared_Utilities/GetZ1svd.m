function [ Z1 ] = GetZ1svd( Y0, Delta )

% function [ Z1 ] = GetZ1svd( Y0, Delta )
% Purpose: Returns Z1 in vectorized form for the Big Problem. Using SVD
%          instead of QR.
% Created:     02.08.2016
% Last change: 12.09.2016

[ n, p ] = size(Delta);

% Do svd of Y0
[ U, ~, V ] = svd( Y0 );
Us = U(:,1:p);
Uperp = U(:,p+1:end);
%Y0 = Us*V';   % if Y0 is not yet orthonormal, this is like doing a projection onto St(n,p)
Y0perp = Uperp;
Q = [ Y0, Y0perp ];

hQ = [ Y0'*Delta,     -Delta'*Y0perp;
      Y0perp'*Delta,    zeros( n-p ) ];

exphQ = expm(hQ);

% Geodesics
Z1 = kron( eye(p,n), Q ) * exphQ(:);

end