function [ Delta_n ] = GetDelta( n, p, Y0, s )

% function [ Delta_n ] = GetDelta( n, p, Y0, s )
% Purpose: Returns a random, unit-norm tangent vector in the tangent space
%          T_{Y0}St(n,p).
% Last change: 18.07.2016

Omega   = rand( s, p, p );
Omega   = Omega - Omega';          % random, p-by-p skew-symmetric matrix
K_0     = rand( s, n, p );         % random, free n-by-p matrix
Delta = Y0*Omega + K_0-Y0*(Y0'*K_0);

% Normalize Delta w.r.t. the "canonical" metric on the tangent space
norm_Delta = GetCanonicalNormDelta( Delta, Y0 );
Delta_n = Delta/norm_Delta;

end