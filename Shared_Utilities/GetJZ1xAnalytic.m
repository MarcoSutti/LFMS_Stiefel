function [ J_Z1_x ] = GetJZ1xAnalytic( E, Q, x )

% function [ J_Z1_x ] = GetJZ1xAnalytic( E, Q, x )
% Purpose: Returns the Jacobian J_Z1_x used in the single shooting method
%          on the Stiefel manifold, calculated analytically.
%          We use the formulation which expresses Z1 as a function of x.
%          x is the column vector containing the free parameters that build
%          up a tangent vector in the tangent space T_{Y0}St(n,p).
%          See notes of 21.09.2016.
%          See also the scripts Derivative_Z1x_wrt_x and
%          Agreement_J_Z1_x_analytic_and_J_Z1_x_FD, in which we
%          numerically verified the correctness of these formulas.
% Created:     22.09.2016
% Last change: 22.09.2016

n = size(Q,1);
p = sqrt(size(E,1));

%--------------------------------------------------------------------------
% Initialize matrices that depend only on n and p
IpZeros = eye( n, p );
P = perfect_shuffle( n-p, p );
% Transformation matrix from blockwise vectorization to ordinary vectorization
H = GetH( n, p );
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

% Get Kronecker representation of the Fréchet derivative of the matrix
% exponential
Khx = GetKA( hx );
sigma_min_K = min(svd(Khx));

% 27.09.2020: Check:
% sigma_min_K - sin(norm(hx,2))/norm(hx,2)

sigma_max_K = max(svd(Khx));

% 27.09.2020: Check:
% sigma_max_K - 1


% Build Jacobian of h(x) wrt x
J_hx = [                        E,        zeros(p^2,p*(n-p));
          zeros(p*(n-p),dimSkewp),              eye(p*(n-p));
          zeros(p*(n-p),dimSkewp),                        -P;
          zeros((n-p)^2,dimSkewp),     zeros((n-p)^2,p*(n-p)) ];

% The Jacobian of Z1 wrt x
J_Z1_x = kron( IpZeros', Q ) * Khx * H * J_hx;

end