n = 2500;
A = getStartMatrix(n, 300);
alpha = 0.5;
%A = delsq(numgrid('S', 32));
%n = size(A, 1);
%b = ones(n, 1);
b = randn(n,1);
b = b / norm(b);
s = 5;
p = 50;
m = s * p + 1;
d = 2 * m;
basisFunc = @mpk;
ctol = 1e-8;

runs = 20;

%define sketch methods
SketchMethods = {...
    'Gaussian',    @Gaussian; ...
    'CountSketch', @CountSketch; ...
    'Rademacher',  @Rademacher ...
    };

%run experiments
iters_all = [];
labels_all = [];
times_all = [];

relErr_runs = {};

for i = 1:size(SketchMethods, 1)
    name = SketchMethods{i, 1};
    genTheta = SketchMethods{i, 2};
    iters = zeros(runs, 1);
    times = zeros(runs, 1);
    Theta0 = genTheta(d, n);

    for j = 1:runs
        Theta = genTheta(d, n); %generate sketch matrix

        if isequal(Theta, Theta0)
            error('Run %d: Theta is SAME as Theta0, stopping program.', j);
        end

        Theta0 = Theta;

        tic;
        [~, beta] = RBGS_GMRES(A, s, p, Theta, basisFunc, b, ctol);
        times(j) = toc;
        iters(j) = length(beta) - 1;

        if strcmp(name, 'Gaussian')
            relErr_runs{end+1} = beta ./ beta(1);
        end

    end

    %collect results
    iters_all = [iters_all; iters];
    times_all = [times_all; times];
    labels_all = [labels_all; repmat({name}, runs, 1)];
end


%BGS_GMRES
tic;
[~, beta] = BGS_GMRES(A, s, p, basisFunc, b, ctol);
time_bgs = toc;
iters_bgs = length(beta) - 1;
%compute relative error
relErr = beta ./ beta(1);



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

% --- RBGS: hollow circles + vertical error bars ---
h = errorbar(0:maxLen-1, relErr_mean, relErr_ci, 'o', ...
    'Color','k','MarkerSize',5,'LineWidth',1.2, ...
    'MarkerFaceColor','w', ...   % hollow circle
    'DisplayName','RBGS-GMRES');
h.LineStyle = 'none';    % no connecting line

% --- BGS: deterministic curve ---
semilogy(0:iters_bgs, relErr, 'r-*', 'LineWidth', 1.5, ...
    'DisplayName','BGS-GMRES');

xlabel('Iteration Number');
ylabel('Relative Residual');
set(gca,'YScale','log');
legend('show','Location','best');
grid on;

exportgraphics(gcf, 'fig/relErr_iter_CI.eps', 'Resolution', 300);




%plot 2: iteration numbers
figure;
boxplot(iters_all, labels_all);
ylabel('Iteration number');
title('Iteration numbers for different sketch methods');
grid on;

ymin = min([iters_all(:); iters_bgs]) - 1;
ymax = max([iters_all(:); iters_bgs]) + 1;
ylim([ymin, ymax]);

%add horizontal line for BGS-GMRES
hold on;
yline(iters_bgs, 'r--', 'LineWidth', 1.5, ...
       'Label', 'BGS_GMRES', 'LabelHorizontalAlignment', 'left');
hold off;
%save plot
exportgraphics(gcf, 'fig/box_iter.eps', 'Resolution', 300);

%plot 3: runtime
figure;
boxplot(times_all, labels_all);
ylabel('Runtime (seconds)');
title('Runtime for different sketch methods');
grid on;
%add horizontal line for BGS_GMRES
hold on;
yline(time_bgs, 'r--', 'LineWidth', 1.5, ...
       'Label', 'BGS_GMRES', 'LabelHorizontalAlignment', 'left');
hold off;
%save plot
exportgraphics(gcf, 'fig/box_runtime.eps', 'Resolution', 300);