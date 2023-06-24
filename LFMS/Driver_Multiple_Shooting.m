%==========================================================================
% Generates Figure 2.5 in my PhD thesis.
% Driver for Multiple Shooting on the Stiefel manifold.
% References: 
% Created:     31.10.2016
% Created:     2016.10.31
% Last change: 2023.03.01

%   Mar 1, 2023:
%       Added save directly into the preprint folder.
%   Nov 22, 2021:
%       Cleanup of comments and other old lines of code.
%   Nov 12, 2021:
%       Added startup file.
%   Nov 18, 2020:
%       Added direct saving into the "Tesi_Marco" folder.
%==========================================================================

% Sets default graphics interpreter, paths and colors.
LFMS_startup;

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Set dimensions of St(n,p)
n = 15;
p = 4;

% Number of subintervals
m = 7;

% Fix stream of random numbers for reproducibility
s = RandStream( 'mt19937ar', 'Seed', 1 );

% % Create random Stiefel matrix Y0
X = eye( n, p );
%X = orth( rand( s, n, p ) );

% Create a random tangent vector Delta in T_{Y0}St(n,p)
distXY = 0.89*pi;
Delta_exact = distXY * GetDelta( n, p, X, s );

global Exact_Length;
Exact_Length = distXY;

% Map the tg vector onto the manifold
[ Y ] = Stiefel_Exp( X, Delta_exact );

%--------------------------------------------------------------------------
% Parameters Multiple Shooting
param.tolSS = 0.5e-16;
param.maxiterSS = 11;

param.verbose = 1;
%--------------------------------------------------------------------------

param.flag  = false;

if p < (n/2)
    % Set up the "baby" formulation:
    [ X_tilde, Y_tilde, Y0perp_tilde, U1 ] = SetBabyProblem( n, p, X, Y );
end

% Define initial guesses:
Sigma_MS_0 = GetStartingGuessSigma( 2*p, p, m, X_tilde, Y_tilde );
    
% Do multiple shooting on the baby problem
[ iterMS, F_MS, norm_F_MS, err_L, Sigma_MS, Delta_rec_baby_MS, param ] = MultipleShootingStiefel_SVD_condensing( 2*p, p, m, Sigma_MS_0, X_tilde, Y_tilde, param );

% We reconstruct the tangent vector at Y0 to St(n,p) from this Delta_rec_baby
A_tilde = X_tilde' * Delta_rec_baby_MS;
B_tilde = Y0perp_tilde' * Delta_rec_baby_MS;
Delta_rec_MS = X*A_tilde + U1*B_tilde;

%--------------------------------------------------------------------------
% Postprocessing of Multiple Shooting
%--------------------------------------------------------------------------
% Convergence plot of multiple shooting:
PlotConvergenceMS( norm_F_MS, err_L );

% All the checks:
MultipleShootingStiefelChecks( 2*p, p, m, Sigma_MS, F_MS, X_tilde, Y_tilde, param.tolSS )
%--------------------------------------------------------------------------
% MS, 18.11.2020:
% Save plot to pdf file in:
fileName_plot = [ 'Plots/Convg_ms_', num2str(n), '_', num2str(p) ];
export_fig(fileName_plot, '-pdf', '-cmyk', '-transparent');
fprintf('--------------------------------------------------------\n');
fprintf('Saved graph to file %s.pdf.\n', fileName_plot);
fprintf('--------------------------------------------------------\n');
%--------------------------------------------------------------------------
% fileName_plot = '../../../../00_Scientific_Research/Shooting_Stiefel_preprint/plots/Convg_ss_15_4';
% export_fig(fileName_plot, '-pdf', '-transparent', '-r300');
% % MS, 2023.03.01:
% % Use option '-dpng' for saving a png image.
% % Use option '-transparent' for transparent background.
% % -r300 is the PPI value, default resolution is low
% fprintf('Saved graph to file %s.pdf.\n', fileName_plot);
% %--------------------------------------------------------------------------

CheckTangentVector( n, p, X, Y, Delta_exact, Delta_rec_MS, param.tolSS )
