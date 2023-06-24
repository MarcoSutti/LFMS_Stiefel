function [ ] = PlotConvergenceSingleShooting( norm_update )

% function [ ] = PlotConvergenceLeapfrog( norm_update )
% Created:     2020.09.09
% Last change: 2023.03.01

%--------------------------------------------------------------------------
% Define some more modern colors:
% yellow = [0.894, 0.671, 0.094];
green = [0.667, 0.706, 0.118];
% gray = [0.325, 0.325, 0.325];
% blue = [0.239, 0.376, 0.655];
% red = [0.827, 0.341, 0.306];
% gray3 = [ 0.255, 0.255, 0.255 ];
% super_light_gray = [ 0.99, 0.99, 0.99 ];
% 
% % Added on 2022.10.07:
% set( 0, 'DefaultAxesXColor', super_light_gray );
% set( 0, 'DefaultAxesYColor', super_light_gray );
% set( 0, 'DefaultAxesZColor', super_light_gray );
% set( 0, 'DefaultTextColor', super_light_gray );

set( 0, 'defaultAxesFontSize', 12 );

%--------------------------------------------------------------------------

% LineWidth of the MarkerEdge:
myMarkerLineWidth = 0.5;

iter = length(norm_update);

% Convergence plot Multiple Shooting
figure

h1 = semilogy( 0:iter-1, norm_update, 'd-', 'Color', green, ...
    'LineWidth', 2, 'MarkerEdgeColor', green, 'MarkerFaceColor', green, ...
    'MarkerSize', 8 );
hold on
h2 = semilogy( 0:0.05:iter-1, (0.27).^(2.^(0:0.05:iter-1) ), '--', ...
    'Color', [0.5 0.5 0.5], 'LineWidth', 1.25 );
grid on

%--------------------------------------------------------------------------
% MS, 2022.10.07:
% Change the color of the grid:
% set( gca,'GridColor', super_light_gray )

% Added on 2022.04.09:
% Questo è il background "esterno" alla figura
% set( gca,'Color', gray3 )
%--------------------------------------------------------------------------

xlabel( 'iteration $k$ of single shooting', 'FontSize', 16 )
% ylabel( '$\|\delta\Delta_{k}\|_{2}$', 'interpreter', 'latex', 'FontSize', 16 )
ylim( [ 0.25*min(norm_update), max(norm_update)*10 ] )
% xlim( [ 0, 8 ] )

drawnow;
h1.MarkerHandle.LineWidth = myMarkerLineWidth;

% Legend
handleLegend = legend( [ h1, h2 ], {'$\|\delta\xi^{(k)}\|_{2}$', 'Quadratic'}, ...
    'FontSize', 14, 'Location', 'NE' );
% handleLegend.EdgeColor = gray3;

drawnow;
lineEntry = findobj(handleLegend.EntryContainer, 'Object', h1 );
entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
entryMarker.LineWidth = myMarkerLineWidth;


% str = sprintf( 'Simple Shooting on St(%d,%d)', n, p );
% title( str, 'interpreter','latex', 'FontSize', 16 )

% dim = [ 0.5, 0, 0, 0.3 ];
% str = sprintf( 'dist($X,Y$)= %0.2f $pi$', distY0Y1/pi );                
% annotation( 'textbox', 'Position', dim, 'String', str, 'FitBoxToText', 'on', ...
%             'BackgroundColor', 'white', ...
%             'interpreter', 'latex', 'FontSize', 16 )
        
end