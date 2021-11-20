function [ Gk ] = GetGkAnalyticSVD( Y0, Delta, n, p, In, Ip, Pnminuspp, Pnn, Pnp, H, IpZeros, ZerosI )

% function [ Gk ] = GetGkAnalyticSVD( Y0, Delta, n, p, In, Ip, Pnminuspp, Pnn, Pnp, H, IpZeros, ZerosI )
% Purpose: Returns the Jacobian Gk used in the multiple shooting method on
%          the Stiefel manifold, calculated using the analytic formulas
%          that we derived for the geodesic formulation which uses the svd
%          (see notes).
%          See also the script Derivative_Z1_and_Z2_wrt_Y0_SVD_12_sept_2016.m
%          in which we numerically verified the correctness of these formulas.
% Created:     12.09.2016
% Last change: 09.11.2016

[ U, S, V ] = svd( Y0 );

% Partition of the svd:
Up    = U(:,1:p);
Uperp = U(:,p+1:end);
Sp    = S(1:p,1:p);
Vp    = V(:,1:p);

Yperp = Uperp;
Q = [ Y0, Yperp ];

hQ = [ IpZeros'*Q'*Delta,   -(ZerosI*Q'*Delta)';
         ZerosI*Q'*Delta,           zeros( n-p ) ];

exphQ = expm(hQ);

KhQ = GetKA( hQ );

F_notE_Y0 = [  kron( Delta', IpZeros' ); 
               kron( Delta',   ZerosI );
               -Pnminuspp * kron( Delta',   ZerosI );
               zeros( (n-p)^2, n^2 ) ];
           
F_notE_Delta = [  kron( Ip, IpZeros'*Q' ); 
                  kron( Ip, ZerosI*Q' );
                  -Pnminuspp*kron( Ip, ZerosI*Q' );
                  zeros( (n-p)^2, n*p ) ];
              
%--------------------------------------------------------------------------
% The Jacobians
J_Uperp_Y0 = -kron( Uperp', Up*inv(Sp)*Vp' ) * Pnp;
J_Q_Y0     = [ eye(n*p); J_Uperp_Y0 ];
J_hQ_Y0    = H * F_notE_Y0 * Pnn * J_Q_Y0;
J_exphQ_Y0 = KhQ * J_hQ_Y0;

J_Z1_Y0 = kron( IpZeros'*exphQ', In ) * J_Q_Y0 + kron( IpZeros', Q ) * J_exphQ_Y0;
J_Z2_Y0 = kron( ( exphQ * hQ * IpZeros )', In ) * J_Q_Y0 ...
          + kron( IpZeros' * hQ', Q ) * J_exphQ_Y0 ...
          + kron( IpZeros', Q * exphQ ) * J_hQ_Y0;
      
J_hQ_Delta    = H * F_notE_Delta;
J_exphQ_Delta = KhQ * J_hQ_Delta;

J_Z1_Delta = kron( IpZeros', Q ) * J_exphQ_Delta;
J_Z2_Delta = kron( IpZeros', Q ) * ( kron( hQ', In ) * J_exphQ_Delta + kron( In, exphQ ) * J_hQ_Delta );

% The block Gk
Gk = [ J_Z1_Y0,  J_Z1_Delta; 
       J_Z2_Y0,  J_Z2_Delta ];

end