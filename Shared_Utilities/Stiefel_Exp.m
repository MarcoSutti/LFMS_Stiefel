function [ Y1 ] = Stiefel_Exp( Y0, Delta )

% function [ Y1 ] = Stiefel_Exp( Y0, Delta )
% Purpose: Compute the Riemannian exponential map on the Stiefel manifold.
% Created:     23.05.2016
% Last change: 13.12.2016
%==========================================================================
% Input arguments:
%   Y0    : base point on St(n,p)
%   Delta : tangent vector in T_{Y0}St(n,p)
% Output arguments:
%   Y1    : Exp^{St}_{Y0}(Delta)
%==========================================================================

[ ~, p ] = size(Y0);

Omega = Y0'*Delta;           % horizontal component of Delta
K = Delta - Y0*Omega;        % normal component of Delta
[ Qe, Re ] = qr( K, 0 );     % "economy-size" qr-decomposition of the normal component of Delta

% Exponential mapping
Y1 = [ Y0, Qe ] * expm( [ Omega, -Re'; Re, zeros(p) ] ) * eye(2*p,p);

end