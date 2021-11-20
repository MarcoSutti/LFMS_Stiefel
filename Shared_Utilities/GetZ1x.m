function [ Z1 ] = GetZ1x( E, Q, x )

% function [ Z1 ] = GetZ1x( E, Q, x )
% Purpose: Returns Z1 in vectorized form for the Big Problem.
%          We use the formulation which expresses Z1 as a function of x.
%          x is the column vector containing the free parameters that build
%          up a tangent vector in the tangent space T_{Y0}St(n,p).
%          See notes of 21.09.2016.
% Created:     22.09.2016
% Last change: 22.09.2016

n = size(Q,1);
p = sqrt(size(E,1));

%--------------------------------------------------------------------------
% Initialize matrices that depend only on n and p
IpZeros = eye( n, p );
%--------------------------------------------------------------------------

% Dimension of Skew(p)
dimSkewp = (1/2)*p*(p-1);

% Extraction of s and d from x:
s = x(1:dimSkewp);    
d = x(dimSkewp+1:end);
    
vecOmega = E*s;

% Matricization of vecOmega
Omega = reshape( vecOmega, [p,p] );

% Matricization of d
D = reshape( d, [n-p,p] );

% Build matrix h(x)
hx = [ Omega,            -D';
           D,    zeros( n-p ) ];

% disp('condition')      
% norm(hx,'fro')^2 < rank(hx)*max(abs(eig(hx)))^2
       
% matexphx
matexphx = expm( hx );

% Geodesics
Z1 = kron( IpZeros', Q ) * matexphx(:);

end