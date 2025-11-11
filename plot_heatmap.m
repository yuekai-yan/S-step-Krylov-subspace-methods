function plot_heatmap(Err, cat, WB_options, AOB_options, ss, matnum)
% Input:
%     Err:   3-D tensor, size = [numWB, numAOB, numS]
%     cat:   category of the error, e.g. relErr, orthErr
%     WB_options, AOB_options:   cell arrays of labels for axes
%     ss:   list of s-values corresponding to the third dimension
%     matnum:   string for dataset name displayed in the title

    [nWB, nAOB, nS] = size(Err);

    % Preprocessing: apply log10 only to positive finite values
    R = log10(Err);

    % Global color limits
    clim = [-16 0];
    % Layout settings: at most 4 subplots per row
    maxCols = 4;
    nCols = min(maxCols, nS);
    nRows = ceil(nS / nCols);

    % Create figure
    figure('Color', 'w');
    tl = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Global title
    sgtitle(sprintf('%s, log10-scale (%s)', cat, matnum), 'FontWeight', 'bold');

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
        vmin_rel = min(Rk(:));
        vmax_rel = max(Rk(:));
        thr_rel  = (vmin_rel + vmax_rel) / 2;
        [m, n] = size(Rk);
        for i = 1:m
          for j = 1:n
            val_rel = Rk(i,j);
            txt_rel = sprintf('%.2f', val_rel);      
            clr = 'w'; if val_rel > thr_rel, clr='k'; end
            text(j,i,txt_rel, 'Color',clr, 'FontWeight','bold',...
                 'HorizontalAlignment','center','VerticalAlignment','middle');
          end
        end
    end

    % Shared colorbar
    cb = colorbar;
    cb.Layout.Tile = 'east';
end
