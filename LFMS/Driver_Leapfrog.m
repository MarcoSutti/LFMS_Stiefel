%==========================================================================
% Driver the leapfrog algorithm described by Noakes, 1998.
% Comparison between Ralf Zimmermann's algorithm for the Stiefel Log and 
% Noakes algorithm.
% Created:     07.06.2021
% Last change: 12.11.2021

%   Nov 12, 2021:
%       Added startup file.
%   Jun 7, 2021:
%       Readded "baby" formulation.
%   Jun 4, 2021:
%       Modified name of file and description.
%   Nov 18, 2020:
%       Added the global variable "Exact_Solution".
%       Commented out the baby problem.
%   July 2, 2020:
%       Added err_L in multiple shooting and its convergence plot.
%   June 29, 2020:
%       Removed the title and the y-axis label from the plots.
%   June 24, 2020:
%       Added 'addpath(genpath('Utilities'))' and 
%       'addpath(genpath('../0_MyUnigeLibrary'))'
%==========================================================================

% Sets default graphics interpreter, paths and colors.
LFMS_startup;

% Global variable to store exact length.
global Exact_Length;
% global Exact_Solution;

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Set dimensions of St(n,p)
n = 12;
p = 3;

% Fix stream of random numbers for reproducibility
s = RandStream( 'mt19937ar', 'Seed', 1 );

% % Create random Stiefel matrix Y0
Y0 = eye( n, p );
% Y0 = orth( rand( s, n, p ) );

% Create a random tangent vector Delta in T_{Y0}St(n,p)
distXY = 0.96*pi;
Delta_exact = distXY * GetDelta( n, p, Y0, s );

Exact_Length = distXY;

% Map the tg vector onto the manifold
[ Y1 ] = Stiefel_Exp( Y0, Delta_exact );

%--------------------------------------------------------------------------
% Parameters Leap-Frog and Multiple Shooting
param.tolLF = 1e-5;
param.maxiterLF = 100;
param.tolSS = 0.5e-14;
param.maxiterSS = 20;

param.verbose = 1;
%--------------------------------------------------------------------------
    
% Ralf Zimmermann's algorithm:
[Delta, k, conv_hist, norm_logV0] = Stiefel_Log_supp(Y0, Y1, 1e-13);
norm(Delta_exact - Delta, 'fro')

[ ~, ~, ~, Delta_rec, param ] = SimpleShootingStiefel_Baby( Y0, Y1, param );

if param.flag==true
    norm(Delta_exact - Delta_rec, 'fro')
    disp('Problem was solved by single shooting.')
    return;
end
param.flag  = false;

if p < (n/2)
    % Set up the "baby" formulation:
    [ Y0_tilde, Y1_tilde, Y0perp_tilde, U1 ] = SetBabyProblem( n, p, Y0, Y1 );
end

% Always start with the minimum number of subintervals
m = 2;

% For the baby problem:
N = 2*p;

while param.flag==false
    
    m = m + 1;
    
    % Initialize F_0
    F_0 = zeros( m*N*p, 1 );

    % Define initial guesses:
    Sigma_LF_0 = GetStartingGuessSigma( N, p, m, Y0_tilde, Y1_tilde );
    
    % 01.07.2020: Compute length at 0th iteration
    [ norm_error_L_omega_0, ~ ] = GetDeltaFromSigma( N, p, m, Sigma_LF_0 );
    
     % Compute F(Sigma):
    %canonical_norm_F = 0;
    for k=1:m-2
        Y0_k    = GetSigmakj( Sigma_LF_0, k, 1, N, p, m );
        Delta_k = GetSigmakj( Sigma_LF_0, k, 2, N, p, m );
        % Calculation of "residuals"
        % NB: the entries of F corresponding to the Stiefel points are
        %     always zero, so we can comment out the following line of code:
        F_0( GetStride( k, 1, N, p, m ) ) = GetZ1svd( Y0_k, Delta_k ) - Sigma_LF_0( GetStride( k+1, 1, N, p, m ) );
        F_0( GetStride( k, 2, N, p, m ) ) = GetZ2svd( Y0_k, Delta_k ) - Sigma_LF_0( GetStride( k+1, 2, N, p, m ) );
    end
    
    norm_F_0 = norm( F_0, 2 ); %canonical_norm_F;

    % Do Leap-Frog on St(2p,p)
    [ iterLF, L_omega, norm_F_LF, norm_error_L_omega, Sigma_LF_k, param ] = LeapFrogStiefel_SVD_2021( N, p, m, Sigma_LF_0, param );
    
end

%%
% 01.07.2020: Concatenate with length at 0th iteration
norm_F_LF = [ norm_F_0, norm_F_LF ];
norm_error_L_omega = [ norm_error_L_omega_0, norm_error_L_omega ];

% Convergence plot of leapfrog:
PlotConvergenceLeapfrog( norm_F_LF, norm_error_L_omega );
