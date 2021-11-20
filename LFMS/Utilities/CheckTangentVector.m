function [ ] = CheckTangentVector( n, p, Y0, Y1, Delta_exact, Delta_reconstructed, tol )

% function [ ] = CheckTangentVector( n, p, Y0, Y1, Delta_exact, Delta_reconstructed, tol )
% Purpose:
% Created:     14.11.2016
% Last change: 26.05.2017

disp('---------------------------------------------------')
disp('               CHECKS TANGENT VECTOR               ')
disp('---------------------------------------------------')

% 1) Is the Delta_reconstructed the same as Delta_exact?
if abs( Delta_reconstructed - Delta_exact ) < tol*ones(n,p)
    disp('Delta_reconstructed is equal to Delta_exact:   OK')
else
    disp('Delta_reconstructed is equal to Delta_exact:   NO :-(')
end

% 2) Is Y1 reached by Delta_reconstructed?
if abs( Stiefel_Exp( Y0, Delta_reconstructed ) - Y1 ) < tol*ones(n,p)
    disp('Y1 is reached by Delta_reconstructed:          OK')
else
    disp('Y1 is reached by Delta_reconstructed:          NO :-(')
end


% 28.10.2016:
CanNorm_Delta_exact = GetCanonicalNormDelta( Delta_exact, Y0 );
CanNorm_Delta_rec = GetCanonicalNormDelta( Delta_reconstructed, Y0 );

% 26.05.2017:
% 3) Is CanonicalNorm_Delta_reconstructed <= CanonicalNorm_Delta_exact?
if (CanNorm_Delta_rec-CanNorm_Delta_exact) <= tol
    disp('CanNorm_Delta_rec <= CanNorm_Delta_exact:      OK')
    formatSpec = '    %5.6f              %5.6f\n';
    fprintf( formatSpec, CanNorm_Delta_rec, CanNorm_Delta_exact )
else
    disp('CanNorm_Delta_rec <= CanNorm_Delta_exact:      NO :-(')
    formatSpec = '    %5.6f              %5.6f\n';
    fprintf( formatSpec, CanNorm_Delta_rec, CanNorm_Delta_exact )
end

end