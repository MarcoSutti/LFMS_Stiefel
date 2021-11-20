function [ iter, L_omega, norm_F, err_L, Sigma_iter, param ] = LeapFrogStiefel_SVD_2021( n, p, m, Sigma_0, param )

% function [ iter, L_omega, norm_F, err_L, Sigma_iter, param ] = LeapFrogStiefel_SVD_2021( n, p, m, Sigma_0, param )
% Purpose: Leapfrog algorithm described by Lyle Noakes, 1996.
% Created:     07.06.2021
% Last change: 07.06.2021

%   Oct 16, 2020:
%       Added Exact_Solution and err_X_k.
%   July 1, 2020:
%       Added Exact_Length as a global variable that stores the exact
%       geodesic distance to be calculated by the leapfrog algorithm.

global Exact_Length;
% global Exact_Solution;

if param.verbose==1
    disp('----------------------------------')
    disp('             LEAP-FROG             ')
    formatSpec = '               m = %1.0d\n';
    fprintf( formatSpec, m );
    disp('----------------------------------')
    fprintf( ' iter     norm_F         err-L\n' );
    disp('----------------------------------')
end

% Initializations:
iter = 0;
L_omega = 0;
norm_F = param.tolLF + 1;
err_L = param.tolLF + 1;
rate  = 1;
F     = zeros( 2*m*n*p, 1 );
Sigma_iter = Sigma_0;
Norm_Deltas = zeros(m-1,1);

while and( norm_F > param.tolLF, iter < param.maxiterLF )   % 03.05.2017: changed stopping criterion, now based on norm_F
    
    iter = iter + 1;
    
    % Cycle over intermediate points
    for k=1:m-2
        %inner_iter = inner_iter + 1;
        Y0_k = GetSigmakj( Sigma_iter, k,   1, n, p, m );
        Y2_k = GetSigmakj( Sigma_iter, k+2, 1, n, p, m );

        % The Midpoint map:
        [ Y1_k, Delta_halved, param ] = MidpointMap( n, p, k, Y0_k, Y2_k, param );
        if param.flag==false
            return;
        end
        
        % Store tangent vector into the array Sigma_iter
        Sigma_iter( GetStride( k, 2, n, p, m ) ) = Delta_halved(:);
        % will get overwritten in next cycle except for k=m-1 (see below)
        %Sigma_iter( GetStride( k+1, 2, n, p, m ) ) = GetZ2svd( Y0_k, Delta_halved );
        
        % Store norm of tangent vector
        Norm_Deltas(k) = GetCanonicalNormDelta( Delta_halved, Y0_k );

        % Store the point in the Sigma array
        Sigma_iter( GetStride( k+1, 1, n, p, m ) ) = Y1_k(:);
        
    end
    
    Delta_second_last = reshape( GetZ2svd( Y0_k, Delta_halved ), [n,p] );
    
    % Store tangent vector into the array Sigma
    Sigma_iter( GetStride( m-1, 2, n, p, m ) ) = Delta_second_last(:);
    
    Y0_second_last = GetSigmakj( Sigma_iter, m-1, 1, n, p, m );
    Norm_Deltas(m-1) = Norm_Deltas(m-2);
    
    % Get last Sigma2:
    Sigma_iter( GetStride( m, 2, n, p, m ) ) = GetZ2svd( Y0_second_last, Delta_second_last );
    
    % Compute F(Sigma):
    %canonical_norm_F = 0;
    for k=1:m-2
        Y0_k    = GetSigmakj( Sigma_iter, k, 1, n, p, m );
        Delta_k = GetSigmakj( Sigma_iter, k, 2, n, p, m );
        % Calculation of "residuals"
        % NB: the entries of F corresponding to the Stiefel points are
        %     always zero, so we can comment out the following line of code:
        F( GetStride( k, 1, n, p, m ) ) = GetZ1svd( Y0_k, Delta_k ) - Sigma_iter( GetStride( k+1, 1, n, p, m ) );
        F( GetStride( k, 2, n, p, m ) ) = GetZ2svd( Y0_k, Delta_k ) - Sigma_iter( GetStride( k+1, 2, n, p, m ) );
        %canonical_norm_F = GetCanonicalNormDelta( reshape( F( GetStride( k, 2, n, p, m ) ), [n,p] ), Y0_k ) + canonical_norm_F;
    end
    
    norm_F(iter) = norm( F, 2 ); %canonical_norm_F;
    
    % Length of piecewise geodesic
    L_omega(iter) = sum(Norm_Deltas);
    
    % Added 02.07.2020:
    err_L(iter) = abs( L_omega(iter) - Exact_Length );
    
    % Added 16.10.2020:
    % Build X_k:
    for k=2:m-1
        X_k((k-2)*n*p+1:(k-1)*n*p) = reshape( GetSigmakj( Sigma_iter, k, 1, n, p, m ), [n*p,1] );
    end
    
    if iter~=1
        rate(iter) = norm_F(iter)/norm_F(iter-1);
    end
    
    if param.verbose==1
        formatSpec = ' %3.2d    %0.4e    %0.4e\n';
        fprintf( formatSpec, iter, norm_F(iter), err_L(iter) )
    end

end

% 03.02.2017: to save the rate of Leapfrog to make the "deterioration of
% rate plot"
param.rateend = rate(end);

end