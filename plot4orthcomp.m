%matnum = 'gen';
%ss = 5:5:40;

[nWB, nAOB, nL] = size(relResAll);
relRes = zeros(nWB, nAOB, nL);
orthLoss = zeros(nWB, nAOB, nL);
idx_ls = zeros(nWB, nAOB, nL);

% AOB_options for RBGS
%AOB_options = ["rCGS2", "rCGS", "RGS", "rMGS"];
%AOB_options = ["rMGS"];
% WB_options for RBGS
%WB_options = ["rCGS", "rCGS2", "RGS", "rMGS", "rWhitening"];
%WB_options = ["rCGS"];

% collect the result for each combination (aob, wb)
mode_case = 2;
switch mode_case
    case 1
        % 1) min(relRes)
        for vs = 1:nL
            s = ss(vs);
            for aob = 1:nAOB
                for wb = 1:nWB
                    relRes_vec = relResAll{wb, aob, vs};
                    orthLoss_vec = orthLossAll{wb, aob, vs};
                    [relRes(wb, aob, vs), idx] = min(relRes_vec);
                    orthLoss(wb, aob, vs) = orthLoss_vec(idx);
                    idx_ls(wb, aob, vs) = (idx) * s;
                end
            end
        end
    case 2
        % 2) < C * min(relRes)
        C = 5;
        for vs = 1:nL
            s = ss(vs);
            for aob = 1:nAOB
                for wb = 1:nWB
                    relRes_vec = relResAll{wb, aob, vs};
                    orthLoss_vec = orthLossAll{wb, aob, vs};
                    minimum = min(relRes_vec(2:end));
                    idx = find(relRes_vec < C * minimum, 1);
                    relRes(wb, aob, vs) =  relRes_vec(idx);
                    orthLoss(wb, aob, vs) = orthLoss_vec(idx);
                    idx_ls(wb, aob, vs) = idx * s;
                end
            end
        end
    case 3
        % 3) use elbow point detection
        for vs = 1:nL
            s = ss(vs);
            for aob = 1:nAOB
                for wb = 1:nWB
                    relRes_vec = relResAll{wb, aob, vs};
                    orthLoss_vec = orthLossAll{wb, aob, vs};
                    y = log10(relRes_vec(:));               
                    valid = isfinite(y);  % remove NaN or Inf entries
                    x = (1:numel(y))';
                    x = x(valid);  y = y(valid);        
                    % line through the first and last points
                    p1 = [x(1), y(1)]; 
                    p2 = [x(end), y(end)];        
                    % perpendicular distance to the endpoint line
                    dist = abs((p2(2)-p1(2))*x - (p2(1)-p1(1))*y + p2(1)*p1(2) - p2(2)*p1(1)) ...
                           / hypot(p2(2)-p1(2), p2(1)-p1(1));                                 
                    % record relRes and orthLoss at the elbow point
                    [~, idx] = max(dist);
                    %[~, idx] = min(relRes_vec);
                    relRes(wb, aob, vs)  = relRes_vec(idx);
                    orthLoss(wb, aob, vs) = orthLoss_vec(idx);
                    idx_ls(wb, aob, vs) = idx * s;
                end
            end
        end
end
% ---------- plot ----------
relRes_log = log10(relRes);
orthLoss_log = log10(orthLoss);
plot_heatmap(relRes_log, "Relative Residual", WB_options, AOB_options, ss, matnum);
plot_heatmap(orthLoss_log, "Loss of Orthogonality", WB_options, AOB_options, ss, matnum);
plot_heatmap(idx_ls, "Index", WB_options, AOB_options, ss, matnum, [min(idx_ls(:)), max(idx_ls(:))]);