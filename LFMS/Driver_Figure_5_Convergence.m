%==========================================================================
% Generates Figure 5 in the leapfrog paper (Figure 3.5 in my PhD thesis).
% Created:     16.10.2020
% Last change: 12.11.2021

%   Nov 12, 2021:
%       Added startup file.
%   Jun 4, 2021:
%       Modified name of file and description.
%   Oct 16, 2020:
%       Created.
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
m_vector = [10,20,50,100];

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
% Parameters Leap-Frog and Shooting
param.tolLF = 1e-15;
param.maxiterLF = 300;
param.tolSS = 0.5e-13;
param.maxiterSS = 20;

param.verbose = 0;
%--------------------------------------------------------------------------

% norm_error_L_omega_matrix = zeros( (m-2)*n*p, 1 )

param.flag  = false;

Legend = cell(length(m_vector),1); 

for iii = 1:length(m_vector)
    
    m = m_vector(iii);
    
    Xstar = zeros( (m-2)*n*p, 1 );
    
    for kkk=1:m-2
        
        Xstar((kkk-1)*n*p+1:kkk*n*p) = reshape( Stiefel_Exp( X, (kkk/(m-1))*Delta_exact ), [n*p,1] );
        
    end
    
    Exact_Solution = Xstar;
    
    F_0 = zeros( 2*m*n*p, 1 );
    
    % Define initial guesses:
    Sigma_LF_0 = GetRandomStartingGuessSigmaLeapfrog( n, p, m, X, Y, param );

%     Sigma_LF_0 = GetStartingGuessSigma( n, p, m, X, Y );
    
    % Build X_k:
    for k=2:m-1
        X_k((k-2)*n*p+1:(k-1)*n*p) = reshape( GetSigmakj( Sigma_LF_0, k, 1, n, p, m ), [n*p,1] );
    end
    
    % Added 16.10.2020:
    err_X_k_0 = norm( Exact_Solution - X_k', 2 );
    
    % 01.07.2020: Compute length at 0th iteration
    [ norm_error_L_omega_0, ~ ] = GetDeltaFromSigma( n, p, m, Sigma_LF_0 );
    
    % Do Leap-Frog on St(2p,p)
    [ iterLF, L_omega, norm_F_LF, norm_error_L_omega, err_X_k, Sigma_LF_k, param ] = LeapFrogStiefel_SVD( n, p, m, Sigma_LF_0, param );
    
    norm_error_L_omega = [ norm_error_L_omega_0, norm_error_L_omega ];
    err_X_k = [ err_X_k_0, err_X_k ];
    
    err_X_k_matrix(:,iii) = err_X_k;
    
    % Save the convergence history for each m:
%     norm_error_L_omega_matrix(:,iii) = norm_error_L_omega';
    
    % This is the Q-convergence definition:
    Q_rate(iii) = err_X_k(end)/err_X_k(end-1);    
    
end
%%
% Plot

% LineWidth of the MarkerEdge:
myMarkerLineWidth = 0.5;

% Stride for controlling the marker and ticks position frequency; i.e.,
% plot a marker and a tick every stride points.
stride = 20;

close all
figure(1)
iter = size(err_X_k,2);

% Convergence plot leapfrog
handle_array(1) = semilogy( 0:iter-1, err_X_k_matrix(:,1), 'd-', 'Color', blue, ...
    'LineWidth', 2, 'MarkerEdgeColor', blue, 'MarkerFaceColor', ...,
    blue, 'MarkerSize', 8, 'MarkerIndices', 1:stride:iter );
hold on
handle_array(2) = semilogy( 0:iter-1, err_X_k_matrix(:,2), 's-', 'Color', yellow, ...
    'LineWidth', 2, 'MarkerEdgeColor', yellow, 'MarkerFaceColor', ...,
    yellow, 'MarkerSize', 8, 'MarkerIndices', 1:stride:iter );
handle_array(3) = semilogy( 0:iter-1, err_X_k_matrix(:,3), 'o-', 'Color', red, ...
    'LineWidth', 2, 'MarkerEdgeColor', red, 'MarkerFaceColor', ...,
    red, 'MarkerSize', 8, 'MarkerIndices', 1:stride:iter );
handle_array(4) = semilogy( 0:iter-1, err_X_k_matrix(:,4), '^-', 'Color', darkgray, ...
    'LineWidth', 2, 'MarkerEdgeColor', darkgray, 'MarkerFaceColor', ...,
    darkgray, 'MarkerSize', 8, 'MarkerIndices', 1:stride:iter );
grid on

xlabel( 'iteration $k$ of leapfrog' )
ylabel( 'err-$k$' )
ylim( [ 1.8e-2, 13 ] )

drawnow;
for i=1:4
    handle_array(i).MarkerHandle.LineWidth = myMarkerLineWidth;
end

% Legend
handleLegend = legend( handle_array, {'10','20','50','100'}, ...
    'FontSize', 14, 'Location', 'SW' );

handleLegend.Title.String = '$m$';

drawnow;
for i=1:4
    lineEntry = findobj(handleLegend.EntryContainer, 'Object', handle_array(i) );
    entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    entryMarker.LineWidth = myMarkerLineWidth;
end

fileName = 'Plots/Convg_LF_12_3_m10_100_maxiter_300';

saveas( gcf, fileName, 'epsc' )
fprintf('-------------------------------------------------------------------------\n');
fprintf('Saved graph to file %s.eps.\n', fileName);
fprintf('-------------------------------------------------------------------------\n');