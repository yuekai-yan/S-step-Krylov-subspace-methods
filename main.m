matnum = 'e20r5000';
%n = 2500;
%A = getStartMatrix(n);
%alpha = 0.5;
%A = delsq(numgrid('S', 32));
%S = load('ML_Geer.mat');
%matrix = S.Problem.A;
matrix = mmread('e20r5000.mtx');

n = size(matrix, 1);
A = @(x) matrix * x;
b = randn(n,1);
b = b / norm(b);
%step size
ss = [5];
m = 1001; 
d = 2 * m;
ctol = 1e-16;
runs = 1;
for j = 1:length(ss)
    s = ss(j);    
    p = (m-1) / s;
    m = s * p + 1;

    %monomial basis
    %basisFunc = @mpk
    
    %Newton basis
    RitzValues = getRitzValues(A, randn(n, 1), s);
    basisFunc = @(Afun, q, s) mpk(Afun, q, s, RitzValues);

    %define sketch methods
    SketchMethods = {...
        'Gaussian',    @Gaussian; ...
        %'CountSketch', @CountSketch; ...
        %'Rademacher',  @Rademacher ...
        };
    
    %run experiments
    iters_all = [];
    labels_all = [];
    times_all = [];
    
    relErr_runs = {};
    orthErr_RBGS = {};
    
    fprintf('Starting RBGS, s = %d \n', s);
    for i = 1:size(SketchMethods, 1)
        name = SketchMethods{i, 1};
        genTheta = SketchMethods{i, 2};
        iters = zeros(runs, 1);
        times = zeros(runs, 1);
        Theta0 = genTheta(d, n);
    
        for j = 1:runs
            Theta = genTheta(d, n); %generate sketch matrix
    
            tic;
            [~, beta, orthErr] = RBGS_GMRES(A, s, p, Theta, basisFunc, b, ctol);
            times(j) = toc;
            iters(j) = length(beta) - 1;
    
            if strcmp(name, 'Gaussian')
                relErr_runs{end+1} = beta;
                orthErr_RBGS{end+1} = orthErr;
            end
    
        end
    
        %collect results
        iters_all = [iters_all; iters];
        times_all = [times_all; times];
        labels_all = [labels_all; repmat({name}, runs, 1)];
    end
    fprintf('RBGS terminated \n');
    
    %BCG_GMRES
    fprintf('Starting BGS, s = %d \n', s);
    tic;
    [~, beta_BCG, orthErr_BGS] = BGS_GMRES(A, s, p, basisFunc, b, ctol);
    time_bgs = toc;
    iters_bgs = length(beta_BCG) - 1;
    %compute relative error
    relErr_BCG = beta_BCG;
    fprintf('BGS terminated \n');
    
    % GMRES matlab
    fprintf('Starting GMRES \n');
    [y,~,~,~,resvec] = gmres(A,b,[],ctol,m);
    fprintf('GMRES terminated \n');
    
    
    %---plot---

    %plot 1: relative error vs. iteration number
    figure; hold on;
    
    % --- RBGS: compute mean and 95% confidence interval ---
    maxLen = max(cellfun(@length, relErr_runs));
    relErr_mat = nan(length(relErr_runs), maxLen);
    for k = 1:length(relErr_runs)
        relErr_mat(k,1:length(relErr_runs{k})) = relErr_runs{k};
    end
    
    relErr_mean = nanmean(relErr_mat,1);
    relErr_std  = nanstd(relErr_mat,[],1);
    nRuns = size(relErr_mat,1);
    
    %95% confidence interval
    relErr_ci = 1.96 * relErr_std / sqrt(nRuns);
    
    % RBGS: hollow circles + vertical error bars
    h = errorbar(0:maxLen-1, relErr_mean, relErr_ci, 'o', ...
        'Color','k','MarkerSize',5,'LineWidth',1, ...
        'MarkerFaceColor','w', ...   % hollow circle
        'DisplayName','RBGS-GMRES');
    h.LineStyle = 'none';    % no connecting line
    
    % BGS
    semilogy(0:iters_bgs, relErr_BCG, 'r-^', 'LineWidth', 1, ...
        'DisplayName','BGS-GMRES');
    
    % GMRES build-in
    semilogy(0:numel(resvec(1:s:end))-1,resvec(1:s:end), 'b-diamond', 'LineWidth', 1, ...
        'DisplayName','GMRES');
    
    xlabel('Iteration Number');
    ylabel('Relative Residual $\| A x - b \|_2 / \| b \|_2$', ...
            'Interpreter', 'latex', 'FontSize', 14);
    set(gca,'YScale','log');
    legend('show','Location','best');
    grid on;
    
    title(sprintf('Relative Residual vs. Iteration Number (s = %d, m = %d, %s)', s, m, matnum));
    
    exportgraphics(gcf, sprintf('fig/%s_relErr_iter_CI_s%d_m%d_reorth_newton.png', matnum, s, m));
    
    
    % Plot 2: orthogonality loss
    figure; hold on;
    
    % --- RBGS: compute mean and 95% confidence interval ---
    orthErr_RBGS_maxLen = max(cellfun(@length, orthErr_RBGS));
    orthErr_RBGS_mat = nan(length(orthErr_RBGS), orthErr_RBGS_maxLen);
    
    for k = 1:length(orthErr_RBGS)
        orthErr_RBGS_mat(k, 1:length(orthErr_RBGS{k})) = orthErr_RBGS{k};
    end
    
    orthErr_RBGS_mean = nanmean(orthErr_RBGS_mat, 1);
    orthErr_RBGS_std  = nanstd(orthErr_RBGS_mat, [], 1);
    orthErr_RBGS_nRuns = size(orthErr_RBGS_mat, 1);
    
    % 95% confidence interval
    orthErr_RBGS_ci = 1.96 * orthErr_RBGS_std ./ sqrt(orthErr_RBGS_nRuns);
    
    % RBGS: hollow circles + vertical error bars
    h = errorbar(0:orthErr_RBGS_maxLen-1, orthErr_RBGS_mean, orthErr_RBGS_ci, 'o', ...
        'Color','k','MarkerSize',5,'LineWidth',1.2, ...
        'MarkerFaceColor','w', ...
        'DisplayName','RBGS-GMRES');
    h.LineStyle = 'none';    % no connecting line
    
    % BGS
    semilogy(0:iters_bgs, orthErr_BGS, 'r-^', 'LineWidth', 1.5, ...
        'DisplayName','BGS-GMRES');
    
    xlabel('Iteration Number');
    ylabel('Orthogonality loss $\|\mathbf{Q}^\top \mathbf{Q} - \mathbf{I}\|_F$', ...
           'Interpreter', 'latex', 'FontSize', 14);
    set(gca,'YScale','log');
    legend('show','Location','best');
    grid on;
    
    title(sprintf('Orthogonality loss vs. Iteration Number (s = %d, m = %d, %s)', s, m, matnum));
    
    
    exportgraphics(gcf, sprintf('fig/%s_orthErr_s%d_m%d_reorth_newton.png', matnum, s, m));
    
    
    
    
    
    
    
        
    
    %plot 3: iteration numbers
    %figure;
    %boxplot(iters_all, labels_all);
    %ylabel('Iteration number');
    %title('Iteration numbers for different sketch methods');
    %grid on;
    
    %ymin = min([iters_all(:); iters_bgs]) - 1;
    %ymax = max([iters_all(:); iters_bgs]) + 1;
    %ylim([ymin, ymax]);
    
    %add horizontal line for BGS-GMRES
    %hold on;
    %yline(iters_bgs, 'r--', 'LineWidth', 1.5, ...
           %'Label', 'BGS_GMRES', 'LabelHorizontalAlignment', 'left');
    %hold off;
    %save plot
    %exportgraphics(gcf, 'fig/box_iter.eps', 'Resolution', 300);
    
    %plot 3: runtime
    %figure;
    %boxplot(times_all, labels_all);
    %ylabel('Runtime (seconds)');
    %title('Runtime for different sketch methods');
    %grid on;
    %add horizontal line for BGS_GMRES
    %hold on;
    %yline(time_bgs, 'r--', 'LineWidth', 1.5, ...
           %'Label', 'BGS_GMRES', 'LabelHorizontalAlignment', 'left');
    %hold off;
    %save plot
    %exportgraphics(gcf, 'fig/box_runtime.eps', 'Resolution', 300);
end