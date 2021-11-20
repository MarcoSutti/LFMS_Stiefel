function [ Y1_k, Delta_halved, param ] = MidpointMap( n, p, k, Y0_k, Y2_k, param )

% function [ Y1_k, Delta_halved, param ] = untitled( n, p, k, Y0_k, Y2_k, param )
% Purpose: Midpoint Map used in the leapfrog algorithm.
% Created:     17.05.2017
% Last change: 31.05.2017

% 1) Get tangent vector using Single Shooting
% 16.01.2017:
[ Y0_tilde, Y2_tilde, Y0perp_tilde, U1 ] = SetBabyProblem( n, p, Y0_k, Y2_k );
Delta_0 = GetStartingGuessDelta( Y0_tilde, Y2_tilde );
[ ~, ~, ~, Delta_Y0Y2_tilde, param ] = SimpleShootingStiefel_BigProblem_Z1x( Y0_tilde, Y2_tilde, Delta_0, param );

% 22.12.2016: If single shooting fails for the current subinterval,
% we exit from the Leapfrog:
if param.flag==false
    if param.verbose==1
        formatSpec = 'Leapfrog: single shooting failed for the subinterval %2.0d. \n Increasing m... \n';
        fprintf( formatSpec, k )
    end
    Delta_halved = 0;
    Y1_k = 0;
    return;
end

% We reconstruct the tangent vector from this solution
A_tilde = Y0_tilde'*Delta_Y0Y2_tilde;
B_tilde = Y0perp_tilde'*Delta_Y0Y2_tilde;
Delta_Y0Y2 = Y0_k*A_tilde + U1*B_tilde;

% 2) Halve the tangent vector
Delta_halved = (1/2) * Delta_Y0Y2;

% 3) Map it on St(n,p)
[ Y1_k ] = Stiefel_Exp( Y0_k, Delta_halved );

end