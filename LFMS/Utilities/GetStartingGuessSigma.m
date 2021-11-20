function [ Sigma ] = GetStartingGuessSigma( n, p, m, Y0, Y1 )

% function [ Sigma ] = GetStartingGuessSigma( n, p, m, Y0, Y1 )
% Purpose: Returns the initial guess for the Multiple Shooting or for the
%          Leap-Frog algorithm.
% Created:     --.--.2016
% Last change: 09.11.2016

% Initializing the solution vector Sigma:
Sigma = zeros( 2*m*n*p, 1 );

% Define starting guesses: loop over the m points
for k=1:m
    if k==1
        Sigma( GetStride( k, 1, n, p, m ) ) = Y0(:);
    else
        Sigma( GetStride( k, 1, n, p, m ) ) = Stiefel_Exp( GetSigmakj( Sigma, k-1, 1, n, p, m ), GetSigmakj( Sigma, k-1, 2, n, p, m ) );
    end
    DeltaFO = 1/(m-k) * ( Y1 - GetSigmakj( Sigma, k, 1, n, p, m ) );
    ProjDeltaFO = ProjTgSpaceStiefel( GetSigmakj( Sigma, k, 1, n, p, m ), DeltaFO );
    Sigma( GetStride( k, 2, n, p, m ) ) = ( norm( DeltaFO )/norm( ProjDeltaFO ) ) * ProjDeltaFO;
end
% Last point, we set it to Y1
Sigma( GetStride( m, 1, n, p, m ) ) = Y1(:);
% Last tangent vector: to get the last Sigma2 we use the derivative of the geodesic, i.e. Z2
Sigma( GetStride( m, 2, n, p, m ) ) = GetZ2svd( GetSigmakj( Sigma, m-1, 1, n, p, m ), GetSigmakj( Sigma, m-1, 2, n, p, m ) );

end


% % Another way
% % Define starting guesses:
% % 1) Sigma11:
% Sigma(1:n*p) = Y0(:);
% 
% % 26.10.2016:
% for k=2:2:2*(m-1)
%     % Intermediate point
%     X = (k/2)/(m-1)*Y1 + (1-(k/2)/(m-1))*Y0;
%     % Projection onto St(n,p)
%     [ U, ~, V ] = svd(X);
%     UVt = U(:,1:p)*V';
%     Sigma(k*n*p+1:(k+1)*n*p) = UVt(:);  
% end
% 
% for k=0:2:2*(m-2)
%     DeltaFO = reshape( Sigma((k+2)*n*p+1:((k+2)+1)*n*p) - Sigma(k*n*p+1:(k+1)*n*p), [n,p] );
%     alpha = norm( DeltaFO );
%     ProjDeltaFO = ProjTgSpaceStiefel( reshape( Sigma(k*n*p+1:(k+1)*n*p), [n,p] ), DeltaFO );
%     beta = norm( ProjDeltaFO );
%     Sigma((k+1)*n*p+1:(k+2)*n*p) = (alpha/beta) * ProjDeltaFO;
% end
% k=2*(m-1);
% % NB: to get the last Sigma2 we use the derivative of the geodesic, i.e. Z2
% Sigma((k+1)*n*p+1:(k+2)*n*p) = GetZ2( reshape( Sigma((k-2)*n*p+1:(k-1)*n*p), [n,p]), reshape( Sigma((k-1)*n*p+1:k*n*p), [n,p]) );
