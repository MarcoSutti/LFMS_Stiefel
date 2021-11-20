function [ Z2 ] = GetZ2svd( Y0, Delta )

% function [ Z2 ] = GetZ2svd( Y0, Delta )
% Created:     02.09.2016
% Last change: 12.09.2016

[ n, p ] = size(Delta);

% Do svd of Y0
[ U, ~, ~ ] = svd( Y0 );
%Us = U(:,1:p); % 15.02.2017
Uperp = U(:,p+1:end);
%Y0 = Us*V';   % if Y0 is not yet orthonormal, this is like doing a projection onto St(n,p)
Y0perp = Uperp;
Q = [ Y0, Y0perp ];

hQ = [ Y0'*Delta,     -Delta'*Y0perp;
      Y0perp'*Delta,    zeros( n-p ) ];

exphQ = expm(hQ);

% Tangent vector:
Z2 = kron( eye(p,n) * hQ', Q ) * exphQ(:);

end