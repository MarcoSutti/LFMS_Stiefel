function [ ] = PlotConvergenceLeapfrog( norm_F, err_L )

% function [ ] = PlotConvergenceLeapfrog( norm_F, err_L )
% Created:     31.05.2017
% Last change: 29.06.2020

%   June 29, 2020:
%       Removed the title and the y-axis label from the plots.

%--------------------------------------------------------------------------
% Define some more modern colors:
% yellow = [0.894, 0.671, 0.094];
green = [0.667, 0.706, 0.118];
% gray = [0.325, 0.325, 0.325];
blue = [0.239, 0.376, 0.655];
% red = [0.827, 0.341, 0.306];
%--------------------------------------------------------------------------

% LineWidth of the MarkerEdge:
myMarkerLineWidth = 0.5;

% Stride for controlling the marker and ticks position frequency; i.e.,
% plot a marker and a tick every stride points.
stride = 5;

iter = length(err_L);

% Convergence plot leapfrog
handle_array(1) = semilogy( 0:iter-1, err_L, '^-', 'Color', green, ...
    'LineWidth', 2, 'MarkerEdgeColor', green, 'MarkerFaceColor', ...,
    green, 'MarkerSize', 8, 'MarkerIndices', 1:stride:iter );
hold on
handle_array(2) = semilogy( 0:iter-1, norm_F, 'd-', 'Color', blue, ...
    'LineWidth', 2, 'MarkerEdgeColor', blue, 'MarkerFaceColor', ...,
    blue, 'MarkerSize', 8, 'MarkerIndices', 1:stride:iter );
handle_array(3) = semilogy( 0:iter-1, 0.03*(1/1.3).^(0:iter-1), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.25 );
grid on
xlabel( 'iteration $k$ of leapfrog', 'FontSize', 16 )
% ylabel( '$| L_{k} - L^{*} |$,~$\| F(\Sigma_{k}) \|_{2}$', 'FontSize', 16 ) %, $\| F(\Sigma) \|_{2}$'
% xlim( [ 1, iter ] )
% ylim( [ min(min(norm_F),min(err_L)), max(max(norm_F),max(err_L))*3 ] )
ylim( [ 1e-6, 10 ] )

% hL = legend( [ h1, h2, h3 ], {'$| L_{k} - L^{*} |$', '$\| F(\Sigma_{k}) \|_{2}$', 'Linear'}, ...
%     'FontSize', 14, 'Location', 'NE' );

drawnow;
for i=1:2
    handle_array(i).MarkerHandle.LineWidth = myMarkerLineWidth;
end

% Legend
handleLegend = legend( handle_array, ...
    {'$| L_{k} - L^{*} |$', ...
    '$\| F(\Sigma_{k}) \|_{2}$', ...
    'Linear'}, ...
    'FontSize', 16, 'Location', 'NE' );

drawnow;
for i=1:2
    lineEntry = findobj(handleLegend.EntryContainer, 'Object', handle_array(i) );
    entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    entryMarker.LineWidth = myMarkerLineWidth;
end

% str = sprintf( 'Leapfrog on St(%d,%d)', n, p );
% title( str, 'FontSize', 16 )

% dim = [ 0.50, 0, 0, 0.90 ];
% str = sprintf( 'dist($X,Y$) = %0.2f $pi$\n $m = %d$', distY0Y1/pi, m );
% annotation( 'textbox', 'Position', dim, 'String', str, 'FitBoxToText', 'on', ...
%     'BackgroundColor', 'white', ...
%     'interpreter', 'latex', 'FontSize', 16 )
end