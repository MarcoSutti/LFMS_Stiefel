function [ iter, FDelta, norm_update, Delta_rec, param ] = SimpleShootingStiefel_Baby( Y0, Y1, param )

% function [ iter, FDelta, norm_update, Delta_rec, param ] = SimpleShootingStiefel_Baby( Y0, Y1, param )
% This version uses the baby problem, where the matrix exponential of a
% 2p-by-2p matrix is taken, instead of a n-by-n matrix.
% The formulation used is the one presented in Rentmeesters thesis, section
% 5.3.
% Created:     20.07.2016
% Last change: 19.01.2017

[ n, p ] = size(Y0);

[ Y0_tilde, Y1_tilde, Y0perp_tilde, U1 ] = SetBabyProblem( n, p, Y0, Y1 );

% 16.01.2017: The generation of an initial guess Delta_0 has been separated
%             from the function SimpleShootingStiefel_BigProblem_Z1x
Delta_0 = GetStartingGuessDelta( Y0_tilde, Y1_tilde );

% 13.10.2016:
[ iter, FDelta, norm_update, Delta_k, param ] = SimpleShootingStiefel_BigProblem_Z1x( Y0_tilde, Y1_tilde, Delta_0, param );

% We reconstruct the tangent vector from this solution
A_tilde = Y0_tilde'*Delta_k;
B_tilde = Y0perp_tilde'*Delta_k;
Delta_rec = Y0*A_tilde + U1*B_tilde;

end