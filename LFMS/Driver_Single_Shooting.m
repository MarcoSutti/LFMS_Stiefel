%==========================================================================
% Generates Figure 2.3 in my PhD thesis.
% Driver for the Simple Shooting Method on the Stiefel Manifold
% using the Kronecker representation of the Frechet derivative of the
% matrix exponential.
% This version uses the baby problem, where the matrix exponential of a
% 2p-by-2p matrix is taken, instead of a n-by-n matrix.
% The formulation used is the one presented in Rentmeesters thesis, section
% 5.3.
%   _   _          _          _   _  _
%  |     |        |            | |    |
%  |  M  |        |  A    -r'  | | Ip |
%  |     | = expm |            | |    |
%  |  N  |        |  r     Op  | | Op |
%  |_   _|        |_          _| |_  _|
   
% Created:     20.07.2016
% Last change: 2022.10.07

%   Oct 7, 2022:
%       Added lines of codes to save a transparent png version of the plot,
%       in order to generate the plot for the talk in the Taipei Postdoc
%       Seminar.
%   Nov 22, 2021:
%       Cleanup of comments and other old lines of code.
%   Nov 12, 2021:
%       Added startup file.
%==========================================================================

% Sets default graphics interpreter, paths and colors.
LFMS_startup;

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Set dimensions of St(n,p)
n = 15;
p = 4;

% Fix stream of random numbers for reproducibility
s = RandStream( 'mt19937ar', 'Seed', 1 );

% Create random Stiefel matrix Y0
% X = rand( n, p );
% [ Y0, ~ ] = qr( X, 0 );   % get the orthogonal factor of X
Y0 = eye( n, p );
%Y0 = orth( rand( n, p ) );
Y0perp = null(Y0');    % The columns of Y0perp span the orthogonal complement to the subspace span(Y0)

% Create a random tangent vector Delta in T_{Y0}St(n,p)
distY0Y1 = 0.75*pi;
Delta_exact = distY0Y1 * GetDelta( n, p, Y0, s );

%load( 'Y0_and_Delta_20_09.mat' );

% Map the tg vector onto the manifold
[ Y1 ] = Stiefel_Exp( Y0, Delta_exact );

param.tolSS = 1e-16;
param.maxiterSS = 9;
param.verbose = 1;
%--------------------------------------------------------------------------

if p < n/2
    [ iter, FDelta, norm_update, Delta_rec, param ] = SimpleShootingStiefel_Baby(  Y0, Y1, param );
else
    Delta_0 = GetStartingGuessDelta( Y0, Y1 );
    [ iter, FDelta, norm_update, Delta_rec, param ] = SimpleShootingStiefel_BigProblem_Z1x( Y0, Y1, Delta_0, param );
end

if param.flag==1
    disp( 'Single shooting worked well.' )
else
    disp( 'Single shooting failed.' )
end

%--------------------------------------------------------------------------
% Postprocessing of Single Shooting
%--------------------------------------------------------------------------
% All the checks:
% SimpleShootingStiefelChecks( Delta_rec, FDelta, Y0, Y1, Delta_exact, param.tolSS )

PlotConvergenceSingleShooting( norm_update )        
% export_fig Plots/Convg_ss_15_4.pdf -pdf -cmyk -transparent;

%--------------------------------------------------------------------------
% SAVE IMAGE TO PNG FILE
% MS, added 2022.10.07.
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save plot                          |\n');
fprintf('+--------------------------------------------------------------+\n');
% MS, 2022.10.07:
% Save plot to png file nella cartella 20221012_Talk_Postdoc_Seminar
% fileName_plot = '../../../../00_NCTS/20221012_Talk_Postdoc_Seminar/plots/Convg_ss_15_4';
fileName_plot = '../../../../00_Scientific_Research/Shooting_Stiefel_preprint/plots/Convg_ss_15_4';
export_fig(fileName_plot, '-pdf', '-transparent', '-r300');
% MS, 2023.03.01:
% Use option '-dpng' for saving a png image.
% Use option '-transparent' for transparent background.
% -r300 is the PPI value, default resolution is low
fprintf('--------------------------------------------------------\n');
fprintf('Saved graph to file %s.eps.\n', fileName_plot);
fprintf('--------------------------------------------------------\n');
%--------------------------------------------------------------------------
