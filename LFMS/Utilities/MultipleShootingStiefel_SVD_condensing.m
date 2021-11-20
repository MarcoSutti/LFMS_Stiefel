function [ iter, F, norm_F, err_L, Sigma, Delta_recovered, param ] = MultipleShootingStiefel_SVD_condensing( n, p, m, Sigma_0, Y0, Y1, param )

% function [ iter, F, norm_F, err_L, Sigma, Delta_recovered, param ] = MultipleShootingStiefel_SVD_condensing( n, p, m, Sigma_0, Y0, Y1, param )
% Purpose: Encapsulates all that is needed to do the multiple shooting.
%          In this version, the global Jacobian DF is not formed.
%          We use the condensing strategy instead, described in Stoer and
%          Bulirsch, p. 519.
% Created:     14.12.2016
% Last change: 03.05.2017

%   July 2, 2020:
%       Added calculation of the geodesic distance.

global Exact_Length;

if param.verbose==1
    disp('--------------------------------------------------------')
    disp('       MULTIPLE SHOOTING WITH CONDENSING STRATEGY       ')
    disp('--------------------------------------------------------')
    fprintf( ' iter     norm_F         err-L      cond(M)     rate \n' );
    disp('--------------------------------------------------------')
end

%--------------------------------------------------------------------------
% Initialize matrices that depend only on n and p
In = eye(n);
Ip = eye(p);
Pnminuspp = perfect_shuffle( n-p, p );
Pnn = perfect_shuffle( n, n );
Pnp = perfect_shuffle( n, p );
H = GetH( n, p );
IpZeros = eye( n, p );
ZerosI = [zeros(n-p,p), eye(n-p)];
% Matrices C and D for the boundary conditions:
C = [ eye(n*p),   zeros(n*p);  zeros(n*p), zeros(n*p) ];
D = [ zeros(n*p), zeros(n*p);    eye(n*p), zeros(n*p) ];
%--------------------------------------------------------------------------

% Initializations for the while loop:
Sigma  = Sigma_0;
update = zeros( 2*m*n*p, 1 );
F      = zeros( 2*m*n*p, 1 );
Gk     = zeros( 2*n*p, 2*n*p, m-1 );
iter   = 0;
norm_update = param.tolSS + 1;
norm_F      = param.tolSS + 1;
rate        = param.tolSS + 1;
param.flag  = true;

while and( norm_update > param.tolSS, iter < param.maxiterSS )
    
    iter = iter + 1;
    
    prodAllGk = 1;
    
    % Cycle over the subintervals
    for k=1:m-1
        Y0_k    = GetSigmakj( Sigma, k, 1, n, p, m );
        Delta_k = GetSigmakj( Sigma, k, 2, n, p, m );
        
        % Added 02.07.2020:
        % Store norm of tangent vector
        Norm_Deltas(k) = GetCanonicalNormDelta( Delta_k, Y0_k );
        
        % Residuals:
        F( GetStride( k, 1, n, p, m ) ) = GetZ1svd( Y0_k, Delta_k ) - Sigma( GetStride( k+1, 1, n, p, m ) );
        F( GetStride( k, 2, n, p, m ) ) = GetZ2svd( Y0_k, Delta_k ) - Sigma( GetStride( k+1, 2, n, p, m ) );
        
        % Compute Jacobians in closed-form expression
        Gk(:,:,k) = GetGkAnalyticSVD( Y0_k, Delta_k, n, p, In, Ip, Pnminuspp, Pnn, Pnp, H, IpZeros, ZerosI );
        prodAllGk = Gk(:,:,k) * prodAllGk;
    end
    
    % Add the BCs to F:
    % BC at the beginning of the interval
    F( GetStride( m, 1, n, p, m ) ) = Sigma( GetStride( 1, 1, n, p, m ) ) - Y0(:);
    % BC at the end of the interval
    F( GetStride( m, 2, n, p, m ) ) = Sigma( GetStride( m, 1, n, p, m ) ) - Y1(:);
    
    % ***** Condensing strategy *****
    M = C + D*prodAllGk;
    
    % 17.01.2017: We put a check on the condition number of M
    if cond(M) > 1e7  % 03.05.2017: increased to 1e4
        if param.verbose==1
            disp( 'MULTIPLE SHOOTING: M HAS BIG CONDITION NUMBER!!!' )
        end
        param.flag = false;
        break;
    end
    
    % Build RHS w:
    SUM = F( [ GetStride( m-1, 1, n, p, m ), GetStride( m-1, 2, n, p, m ) ] );
    prodGk = 1;
    for k=m-1:-1:2
        prodGk = prodGk * Gk(:,:,k);
        SUM = SUM + prodGk * F( [ GetStride( k-1, 1, n, p, m ), GetStride( k-1, 2, n, p, m ) ] );
    end
    
    w = -( F( [ GetStride( m, 1, n, p, m ), GetStride( m, 2, n, p, m ) ] ) + D * SUM );
    
    % Solve for update1 with LU factorization
    update( [ GetStride( 1, 1, n, p, m ), GetStride( 1, 2, n, p, m ) ] ) = linsolve( M, w );
    
    % Get all the other update's ("waterfallwise")
    for k=2:m
        update( [ GetStride( k, 1, n, p, m ), GetStride( k, 2, n, p, m ) ] ) = Gk(:,:,k-1)*update( [ GetStride( k-1, 1, n, p, m ), GetStride( k-1, 2, n, p, m ) ] ) + F( [ GetStride( k-1, 1, n, p, m ), GetStride( k-1, 2, n, p, m ) ] );
    end
    % ***** end of condensing strategy *****
    
    % Update
    Sigma = Sigma + update;
    
    
    % Added 02.07.2020:
    % Length of piecewise geodesic
    L_omega(iter) = sum(Norm_Deltas);
    
    
    % Monitor the norm of the update and of F
    norm_update(iter) = norm( update, 2 );
    norm_F(iter)      = norm( F, 2 );
    if iter~=1
        rate(iter) = norm_F(iter)/norm_F(iter-1);
    end
    
    % 01.07.2020: Added Exact_Length
    err_L(iter) = abs( L_omega(iter) - Exact_Length );
    
    if param.verbose==1
        if iter==1 % at the first iteration, rate has no sense
            formatSpec = ' %3.2d    %.4e    %.4e    %.3f \n';
            fprintf( formatSpec, iter-1, norm_F(iter), err_L(iter), cond(M) )
        else
            formatSpec = ' %3.2d    %.4e    %.4e    %.3f      %.3f \n';
            fprintf( formatSpec, iter-1, norm_F(iter), err_L(iter), cond(M), rate(iter) )
        end
    end
    
    %----------------------- DEBUGGING -----------------------
    % PlotPointsAndVectors( n, p, m, Y0, Y1, Sigma )
    %----------------------- DEBUGGING -----------------------
end

if iter==param.maxiterSS
    % 17.01.2017: if we reach the maximum number of iterations allowed for
    % Multiple Shooting, we set the flag to 'false'
    if param.verbose==1
        disp( 'MULTIPLE SHOOTING: REACHED MAXIMUM NUMBER OF ITERATIONS!!!' )
    end
    param.flag = false;
    return;
end

% We build the tangent vector at Y0 to St(n,p) (Last change: 15.02.2017)
[ ~, Delta_recovered ] = GetDeltaFromSigma( n, p, m, Sigma );

end