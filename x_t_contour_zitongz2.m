function hC = X_T_contour(ax, X, T, Z, x_ind, t_ind, z_lim, title_str, ...
    edge_color, c_scale, cmap)
if strcmpi(c_scale, 'log'); Z_plot = log10(Z);
else                        Z_plot = Z;         end

hC = contourf(ax, X(t_ind,x_ind), T(t_ind,x_ind), Z_plot(t_ind,x_ind), ...
    'edgecolor', edge_color);

set(ax, 'FontSize', 14, 'FontName', 'Times New Roman');

switch cmap
    case 'black_body'
        color_map = xlsread('black-body-table-float-1024.csv');
        color_map = color_map(:,2:4);           colormap(ax, color_map);
    case 'cool_warm'
        color_map = xlsread('cool-warm-table.csv');
        color_map = color_map(:,2:4)/255;       colormap(ax, color_map);
    case 'green_seq'
        color_map = xlsread('green_seq.csv');   color_map = color_map/255;
        colormap(ax, color_map);
    case 'purple_seq'
        color_map = xlsread('purple_seq.csv');  color_map = color_map/255;
        colormap(ax, color_map);
    case 'warm'
        color_map = xlsread('cool-warm-table.csv');
        color_map = color_map(ceil(0.5*size(color_map,1)):end,2:4)/255;
        colormap(ax, color_map);
    case 'cool'
        color_map = xlsread('cool-warm-table.csv');
        color_map = color_map(floor(0.5*size(color_map,1)):-1:1,2:4)/255;
        colormap(ax, color_map);
    otherwise
        try
            colormap(ax, cmap);
        catch
            colormap(ax, 'gray');
        end
end

if ~isempty(z_lim)
    if strcmpi(c_scale, 'log'); caxis(ax, log10(z_lim));
    else                        caxis(ax, z_lim);           end
end

h_cb = colorbar(ax);
if strcmpi(c_scale, 'log')
    h_cb.TickLabels = num2cell(10.^h_cb.Ticks);
end

xlabel(ax, '$ X\ (m) $', 'interpreter', 'latex'); 
ylabel(ax, '$ t\ (s) $', 'interpreter', 'latex'); 
title(ax, title_str, 'FontSize', 16, 'interpreter', 'latex');