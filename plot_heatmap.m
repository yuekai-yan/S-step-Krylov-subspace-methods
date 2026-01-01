function plot_heatmap(R, cat, WB_options, AOB_options, ss, matnum, clim)
% Input:
%     Err:   3-D tensor, size = [numWB, numAOB, numS]
%     cat:   category of the error, e.g. relErr, orthErr
%     WB_options, AOB_options:   cell arrays of labels for axes
%     ss:   list of s-values corresponding to the third dimension
%     matnum:   string for dataset name displayed in the title

    if nargin < 7
        clim = [-16 0];
    end

    [nWB, nAOB, nS] = size(R);
    % Layout settings: at most 4 subplots per row
    maxCols = 4;
    nCols = min(maxCols, nS);
    nRows = ceil(nS / nCols);

    % Create figure
    figure('Color', 'w');
    tl = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Global title
    sgtitle(sprintf('%s (%s)', cat, matnum), 'FontWeight', 'bold');

    for k = 1:nS
        nexttile;

        % Current slice
        Rk = R(:,:,k);

        % Heatmap (use AlphaData to hide invalid values)
        hImg = imagesc(Rk, clim);
        colormap(parula);

        % subplot
        set(gca, 'XTick', 1:nAOB, 'XTickLabel', AOB_options, ...
                 'YTick', 1:nWB, 'YTickLabel', WB_options, ...
                 'TickDir', 'out', 'Layer', 'top', 'FontName', 'Helvetica');
        xtickangle(35);
        axis tight;
        box on;
        title(sprintf('s = %d', ss(k)));

        % Add numeric values inside cells
        [m, n] = size(Rk);
        for i = 1:m
          for j = 1:n
            val_rel = Rk(i,j);
            if strcmp(cat, "Index")
                txt_rel = sprintf('%.0f', val_rel);
            else
                txt_rel = sprintf('%.2f', val_rel);
            end
            clr = 'k'; 
            text(j,i,txt_rel, 'Color',clr, 'FontWeight','bold',...
                 'HorizontalAlignment','center','VerticalAlignment','middle');
          end
        end
    end

    % Shared colorbar
    cb = colorbar;
    cb.Layout.Tile = 'east';
end
