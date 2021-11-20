%==========================================================================
% Generates Figures 4.9 and 4.10 in the thesis.
% Driver for Multiple Shooting which uses an initial guess provided by the
% leapfrog algorithm described by Noakes, 1998.
% References: 
% Created:     31.10.2016
% Last change: 12.11.2021

%   Nov 12, 2021:
%       Added startup file.
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

% Global variables to store exact length and exact solution:
global Exact_Length;
global Exact_Solution;

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Set dimensions of St(n,p)
n = 12;
p = 3;

% Fix stream of random numbers for reproducibility
s = RandStream( 'mt19937ar', 'Seed', 1 );

% % Create random Stiefel matrix Y0
X = eye( n, p );
%X = orth( rand( s, n, p ) );

% Create a random tangent vector Delta in T_{Y0}St(n,p)
distXY = 0.95*pi;
Delta_exact = distXY * GetDelta( n, p, X, s );

Exact_Length = distXY;

% Map the tg vector onto the manifold
[ Y ] = Stiefel_Exp( X, Delta_exact );

%--------------------------------------------------------------------------
% Parameters Leap-Frog and Multiple Shooting
param.tolLF = 1e-3;
param.maxiterLF = 100;
param.tolSS = 0.5e-14;
param.maxiterSS = 20;

param.verbose = 1;
%--------------------------------------------------------------------------
    
% [ ~, ~, ~, Delta_rec, param ] = SimpleShootingStiefel_Baby( Y0, Y1, param );
%
% if param.flag==true
%     disp('Problem was solved by single shooting.')
%     return;
% end
param.flag  = false;

% if p < (n/2)
%     % Set up the "baby" formulation:
%     [ X_tilde, Y_tilde, Y0perp_tilde, U1 ] = SetBabyProblem( n, p, X, Y );
% end

% Always start with the minimum number of subintervals
m = 3;

F_0 = zeros( 2*m*n*p, 1 );

% For the baby problem:
N = n;

while param.flag==false
    
    m = m + 1;
    
    % MS, 18.11.2020: Define exact solution:
    Xstar = zeros( (m-2)*N*p, 1 );
    
    for kkk=1:m-2

        Xstar((kkk-1)*N*p+1:kkk*N*p) = reshape( Stiefel_Exp( X, (kkk/(m-1))*Delta_exact ), [N*p,1] );
        
    end
    
    Exact_Solution = Xstar;
    
    % Define initial guesses:
    Sigma_LF_0 = GetStartingGuessSigma( N, p, m, X, Y );

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
    [ iterLF, L_omega, norm_F_LF, norm_error_L_omega, err_X_k, Sigma_LF_k, param ] = LeapFrogStiefel_SVD( N, p, m, Sigma_LF_0, param );
    
end

%%
% 01.07.2020: Concatenate with length at 0th iteration
norm_F_LF = [ norm_F_0, norm_F_LF ];
norm_error_L_omega = [ norm_error_L_omega_0, norm_error_L_omega ];

% Convergence plot of leapfrog:
PlotConvergenceLeapfrog( norm_F_LF, norm_error_L_omega );
% export_fig Plots/Convg_lf_12_3.pdf -pdf -cmyk -transparent;
%--------------------------------------------------------------------------
% MS, 18.11.2020:
% Save plot to png file nella cartella Tesi_Marco
% fileName = [ '../../../Dropbox/0_UNIGE/Tesi_Marco/images/simulations/Convg_lf_', num2str(n), '_', num2str(p) ];
fileName = [ 'Plots/Convg_lf_', num2str(n), '_', num2str(p) ];
export_fig(fileName, '-pdf', '-cmyk', '-transparent');
fprintf('Saved graph to file %s.pdf.\n', fileName);
%--------------------------------------------------------------------------

% %--------------------------------------------------------------------------
% % Save the plot to file
% fileName = 'Leapfrog_Plots/Convg_lf_12_3';
% pause(0.5)
% saveas( gcf, fileName, 'epsc' );
% fprintf('+-------------------------------------------------+\n');
% fprintf('Saved graph to file %s.eps.\n', fileName);


% Do Multiple Shooting using the Sigma's obtained with the leapfrog algorithm
Sigma_MS_0 = Sigma_LF_k;


% % 24.05.2017: do the alternative shooting with the broken geodesic
% [ iterBGS, F_BGS, norm_F_BGS, norm_update, Sigma_BGS, Delta_rec_baby_BGS, param ] = BrokenGeodesicShooting( N, p, m, Sigma_MS_0, Y1_tilde, param );
%
% % We reconstruct the tangent vector at Y0 to St(n,p) from this Delta_rec_baby
% A_tilde = Y0_tilde' * Delta_rec_baby_BGS;
% B_tilde = Y0perp_tilde' * Delta_rec_baby_BGS;
% Delta_rec_BGS = Y0*A_tilde + U1*B_tilde;
%
% % Convergence plot of broken geodesic shooting:
% %PlotConvergenceBGS( iterBGS, norm_F_BGS, norm_update, distY0Y1, n, p, m );

%% Do multiple shooting on the baby problem
[ iterMS, F_MS, norm_F_MS, err_L, Sigma_MS, Delta_rec_MS, param ] = MultipleShootingStiefel_SVD_condensing( N, p, m, Sigma_MS_0, X, Y, param );

% % We reconstruct the tangent vector at Y0 to St(n,p) from this Delta_rec_baby
% A_tilde = X_tilde' * Delta_rec_baby_MS;
% B_tilde = Y0perp_tilde' * Delta_rec_baby_MS;
% Delta_rec_MS = X*A_tilde + U1*B_tilde;

%--------------------------------------------------------------------------
% Postprocessing of Multiple Shooting
%--------------------------------------------------------------------------
% Convergence plot of multiple shooting:
PlotConvergenceMS( norm_F_MS, err_L );
% MS, 18.11.2020:
% Save plot to png file nella cartella Tesi_Marco
% fileName = [ '../../../Dropbox/0_UNIGE/Tesi_Marco/images/simulations/Convg_ms_', num2str(n), '_', num2str(p) ];
fileName = [ 'Plots/Convg_ms_', num2str(n), '_', num2str(p) ];
export_fig(fileName, '-pdf', '-cmyk', '-transparent');
fprintf('Saved graph to file %s.pdf.\n', fileName);
%--------------------------------------------------------------------------
% All the checks:
MultipleShootingStiefelChecks( N, p, m, Sigma_MS, F_MS, X, Y, param.tolSS )
%--------------------------------------------------------------------------


% %--------------------------------------------------------------------------
% % Save the plot to file
% fileName = 'Leapfrog_Plots/Convg_ms_12_3';
% pause(0.5)
% saveas( gcf, fileName, 'epsc' );
% fprintf('+-------------------------------------------------+\n');
% fprintf('Saved graph to file %s.eps.\n', fileName);


%CheckTangentVector( n, p, X, Y, Delta_exact, Delta_rec_BGS, param.tolSS )
CheckTangentVector( n, p, X, Y, Delta_exact, Delta_rec_MS, param.tolSS )


% --------------------------------------------------------------------------
% % 10.05.2017: plot of norm_F also with intermediate iterations:
% figure(1)
% % semilogy( 1:length(norm_F_LF), norm_F_LF, 'ko-' )
% % semilogy( 1:length(norm_F_LF_RB), norm_F_LF_RB, 'rd-' )
% semilogy( 1:length(norm_error_L_omega), norm_error_L_omega, 'bo-' )
% hold on
% %semilogy( 1:length(norm_error_L_omega_RB), norm_error_L_omega_RB, 'rd-' )
% xlabel( 'inner iterations', 'interpreter', 'latex', 'FontSize', 16 )
% ylabel( '$| L(\omega_{a}) - L^{*} |$', 'interpreter', 'latex', 'FontSize', 16 )
% str = sprintf( 'Leap-Frog on St(%d,%d)', n, p );
% title( str, 'interpreter','latex', 'FontSize', 16 )
%
% dim = [ 0.5, 0, 0, 0.3 ];
% str = sprintf( 'dist($Y_{0},Y_{1}$)= %0.2f $pi$\n $m = %d$', distY0Y1/pi, m );
% annotation( 'textbox', 'Position', dim, 'String', str, 'FitBoxToText', 'on', ...
%             'BackgroundColor', 'white', ...
%             'interpreter', 'latex', 'FontSize', 16 )
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% figure(2)
% semilogy( 1:length(norm_F_LF), norm_F_LF, 'bo-' )
% hold on
% %semilogy( 1:length(norm_F_LF_RB), norm_F_LF_RB, 'rd-' )
% xlabel( 'inner iterations', 'interpreter', 'latex', 'FontSize', 16 )
% ylabel( '$\| F(\Sigma) \|_{2}$', 'interpreter', 'latex', 'FontSize', 16 )
% str = sprintf( 'Leap-Frog on St(%d,%d)', n, p );
% title( str, 'interpreter','latex', 'FontSize', 16 )
%
% dim = [ 0.5, 0, 0, 0.3 ];
% str = sprintf( 'dist($Y_{0},Y_{1}$)= %0.2f $pi$\n $m = %d$', distY0Y1/pi, m );
% annotation( 'textbox', 'Position', dim, 'String', str, 'FitBoxToText', 'on', ...
%             'BackgroundColor', 'white', ...
%             'interpreter', 'latex', 'FontSize', 16 )
% %--------------------------------------------------------------------------
%sum(param.time_inverting_J)
%sum(param.time_condensing)