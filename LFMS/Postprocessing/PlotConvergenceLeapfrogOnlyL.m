function [ ] = PlotConvergenceLeapfrogOnlyL( norm_error_L_omega_matrix, m_vector )

% function [ ] = PlotConvergenceLeapfrogOnlyL( norm_error_L_omega_matrix )
% Created:     16.10.2020
% Last change: 16.10.2020

%   Oct 16, 2020:
%       Created.
%--------------------------------------------------------------------------

% LineWidth of the MarkerEdge:
% myMarkerLineWidth = 0.5;

% Stride for controlling the marker and ticks position frequency; i.e.,
% plot a marker and a tick every stride points.
stride = 5;

iter = size(norm_error_L_omega_matrix,1);

% Convergence plot leapfrog

semilogy( 0:iter-1, norm_error_L_omega_matrix, 'd-', 'LineWidth', 2, ...
    'MarkerIndices', 1:stride:iter );

grid on
xlabel( 'iteration $k$ of leapfrog', 'FontSize', 16 )
ylabel( '$| L_{k} - L^{*} |$', 'FontSize', 16 )
% ylim( [ 1e-6, 10 ] )

% lgd = legend({num2str(m_vector(1)), num2str(m_vector(2)), num2str(m_vector(3)), num2str(m_vector(4)), ...
%     num2str(m_vector(5)), num2str(m_vector(6)), num2str(m_vector(7))},'Location','NE');
% lgd.Title.String = '$m$';

% 
% drawnow;
% handle_array(1).MarkerHandle.LineWidth = myMarkerLineWidth;
% 
% % Legend
% handleLegend = legend( handle_array, ...
%     {'$| L_{k} - L^{*} |$', ...
%     'Linear'}, ...
%     'FontSize', 14, 'Location', 'NE' );
% 
% drawnow;
% lineEntry = findobj(handleLegend.EntryContainer, 'Object', handle_array(1) );
% entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
% entryMarker.LineWidth = myMarkerLineWidth;

end