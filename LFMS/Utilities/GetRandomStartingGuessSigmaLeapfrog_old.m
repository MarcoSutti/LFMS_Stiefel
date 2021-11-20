function [ Sigma ] = GetRandomStartingGuessSigmaLeapfrog_old( n, p, m, Y0, Y1, param )

% function [ Sigma ] = GetRandomStartingGuessSigmaLeapfrog( n, p, m, Y0, Y1, param )
% Purpose: Returns a random initial guess for the Leap-Frog algorithm.
% Created:     19.10.2020
% Last change: 19.10.2020

% Initializing the solution vector Sigma:
Sigma = zeros( 2*m*n*p, 1 );

Sigma( GetStride( 1, 1, n, p, m ) ) = Y0(:);
Sigma( GetStride( m, 1, n, p, m ) ) = Y1(:);

for k=2:m-1
    % Random points
    Sigma( GetStride( k, 1, n, p, m ) ) = orth( rand( n, p ) );
end

% Find the tangent vectors:
for k=1:m-1
    Y0_k = reshape( Sigma( GetStride( k, 1, n, p, m ) ), [n,p] );
    Y1_k = reshape( Sigma( GetStride( k+1, 1, n, p, m ) ), [n,p] );
    Delta_0_k = GetStartingGuessDelta( Y0_k, Y1_k );
    [ ~, ~, ~, Delta_k, param ] = SimpleShootingStiefel_BigProblem_Z1x( Y0_k, Y1_k, Delta_0_k, param );
    Sigma( GetStride( k, 2, n, p, m ) ) = Delta_k(:);
end

% Last tangent vector: to get the last Sigma2 we use the derivative of the geodesic, i.e. Z2
Sigma( GetStride( m, 2, n, p, m ) ) = GetZ2svd( GetSigmakj( Sigma, m-1, 1, n, p, m ), GetSigmakj( Sigma, m-1, 2, n, p, m ) );

end