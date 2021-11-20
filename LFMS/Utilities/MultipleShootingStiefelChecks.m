function [ ] = MultipleShootingStiefelChecks( n, p, m, Sigma, F, Y0, Y1, tol )

% function [ ] = MultipleShootingStiefelChecks( n, p, m, Sigma, F, Y0, Y1, Delta_exact, tol )
% Purpose: Contains all the checks to be performed on the solution Sigma
%          found by the Multiple Shooting on the Stiefel manifold.
% Created:     01.09.2016
% Last change: 14.11.2016

disp('---------------------------------------------------')
disp('             CHECKS MULTIPLE SHOOTING              ')
disp('---------------------------------------------------')

% 1) First point must coincide with Y0 and last must coincide with Y1
if abs( Sigma( GetStride( 1, 1, n, p, m ) ) - Y0(:) ) < tol*ones(n*p,1)
    disp('Sigma11 is Y0:             OK')
else
    disp('Sigma11 is Y0:             NO :-(')
end
if abs( Sigma( GetStride( m, 1, n, p, m ) ) - Y1(:) ) < tol*ones(n*p,1)
    formatSpec = 'Sigma%0.0d1 is Y1:             OK \n';
    fprintf( formatSpec, m )
else
    formatSpec = 'Sigma%0.0d1 is Y1:             NO :-( \n';
    fprintf( formatSpec, m )
end

% 2) The intermediate points have to be orthonormal:
for k=2:m-1
    if abs( GetSigmakj( Sigma, k, 1, n, p, m )'*GetSigmakj( Sigma, k, 1, n, p, m ) - eye(p) ) < tol*ones(p)
        formatSpec = 'Sigma%0.0d1 in St(n,p):        OK \n';
        fprintf( formatSpec, k )
    else
        formatSpec = 'Sigma%0.0d1 in St(n,p):        NO :-( \n';
        fprintf( formatSpec, k )
    end
end

% 4) Check that F is identically zero (component-wise)
if abs( F ) < tol*ones(2*m*n*p,1)
    disp('F is numerically zero:     OK')
else
    disp('F is numerically zero:     NO :-(')
end

end

