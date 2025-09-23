A = delsq(numgrid('S', 32));
n = size(A, 1);
b = ones(n, 1);
d = 500;
s = 5;
p = 30;
m = s * p + 1;
basisFunc = @mpk;
ctol = 1e-6;

runs = 20;

%define sketch methods
SketchMethods = {...
    'Gaussian',    @Gaussian; ...
    'CountSketch', @CountSketch ...
    };

%run experiments
iters_all = [];
labels_all = [];
times_all = [];

for i = 1:size(SketchMethods, 1)
    name = SketchMethods{i, 1};
    genTheta = SketchMethods{i, 2};
    iters = zeros(runs, 1);
    times = zeros(runs, 1);

    for j = 1:runs
        Theta = genTheta(d, n); %generate sketch matrix
        tic;
        [~, beta] = RBGS_GMRES(A, s, p, Theta, basisFunc, b, ctol);
        times(j) = toc;
        iters(j) = length(beta) - 1;
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


%plot 1: iteration numbers
figure;
boxplot(iters_all, labels_all);
ylabel('Iteration number');
title('Iteration numbers for different sketch methods');
grid on;
%add horizontal line for BGS-GMRES
hold on;
yline(iters_bgs, 'r--', 'LineWidth', 1.5, ...
       'Label', 'BGS_GMRES', 'LabelHorizontalAlignment', 'left');
hold off;
%save plot
exportgraphics(gcf, 'fig/box_iter.eps', 'Resolution', 300);


%plot 2: runtime
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