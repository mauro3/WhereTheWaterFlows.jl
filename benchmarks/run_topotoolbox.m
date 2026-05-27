clear; clc;

%We compare the time required to compute single-flow-direction routing, flow accumulation, and drainage basin labels from the same DEM. 
% For WWF this is one waterflows call; for TopoToolbox this is FLOWobj + flowacc + drainagebasins.

threads_env = getenv('TOPO_THREADS');
if isempty(threads_env)
    nthreads_req = 1;
else
    nthreads_req = str2double(threads_env);
    if isnan(nthreads_req) || nthreads_req < 1
        error('TOPO_THREADS must be a positive integer')
    end
end

dataset_env = getenv('TOPO_DATASET');
if isempty(dataset_env)
    dataset = 'large';
else
    dataset = lower(strtrim(dataset_env));
end
if ~strcmp(dataset, 'small') && ~strcmp(dataset, 'large')
    error('TOPO_DATASET must be "small" or "large"')
end

runs_env = getenv('TOPO_RUNS');
if isempty(runs_env)
    nruns = 6;
else
    nruns = str2double(runs_env);
    if isnan(nruns) || nruns < 1 || mod(nruns, 1) ~= 0
        error('TOPO_RUNS must be a positive integer')
    end
end

maxNumCompThreads(nthreads_req);
nthreads = maxNumCompThreads;

fprintf('MATLAB threads: %d\n', nthreads);
fprintf('Dataset: %s\n', dataset);
fprintf('Runs: %d\n', nruns);

% Reproducible paths
root = fileparts(mfilename('fullpath'));

topodir = fullfile(root, 'topotoolbox3');

if ~exist(topodir, 'dir')
    error('TopoToolbox not found. Clone it into benchmarks/: git clone https://github.com/TopoToolbox/topotoolbox3.git topotoolbox3')
end

addpath(genpath(topodir))

datadir = fullfile(root, 'data_raw');
outdir  = fullfile(root, 'outputs', 'topotoolbox');

if ~exist(datadir, 'dir')
    mkdir(datadir);
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

demfile_small = fullfile(datadir, 'swissalti3d_tile.tif');
demfile_large = fullfile(datadir, 'swissalti3d_tilelarge.tif');

if strcmp(dataset, 'large') && ~exist(demfile_large, 'file')
    error(['Large DEM not found: ' demfile_large newline ...
           'Create it first by running WWF once in real-large mode, e.g.:' newline ...
           'julia --project run_wwf.jl --mode real --dataset large --runs 1']);
end

if strcmp(dataset, 'small')
    demfile = demfile_small;
else
    demfile = demfile_large;
end

% -------------------------------------------------------
% Load DEM
% -------------------------------------------------------

DEM = GRIDobj(demfile);

npixels = nnz(~isnan(DEM.Z));


% -------------------------------------------------------
% Run TopoToolbox benchmark
% -------------------------------------------------------

runtime_flowdir_s = zeros(nruns, 1);
runtime_flowacc_s = zeros(nruns, 1);
runtime_basins_s  = zeros(nruns, 1);
runtime_total_s   = zeros(nruns, 1);

for i = 1:nruns
    fprintf('TopoToolbox run %d / %d\n', i, nruns);

    tic
    FD = FLOWobj(DEM);
    runtime_flowdir_s(i) = toc;

    tic
    A = flowacc(FD);
    runtime_flowacc_s(i) = toc;

    tic
    DB = drainagebasins(FD);
    runtime_basins_s(i) = toc;

    runtime_total_s(i) = ...
        runtime_flowdir_s(i) + ...
        runtime_flowacc_s(i) + ...
        runtime_basins_s(i);

    fprintf('  runtime: %.3f s\n', runtime_total_s(i));
end

features = "FLOWobj + flowacc + drainagebasins";

% -------------------------------------------------------
% Save benchmark CSV
% -------------------------------------------------------

if strcmp(dataset, 'large')
    warmup_drop = min(1, nruns - 1);
else
    warmup_drop = min(2, nruns - 1);
end
mean_runtime = mean(runtime_total_s((warmup_drop + 1):end));
std_runtime  = std(runtime_total_s((warmup_drop + 1):end));

Tnew = table( ...
    "topotoolbox", ...
    features, ...
    npixels, ...
    nthreads, ...
    nruns, ...
    mean_runtime, ...
    std_runtime, ...
    'VariableNames', { ...
        'method', ...
        'features', ...
        'npixels', ...
        'nthreads', ...
        'nruns', ...
        'runtime_mean_s', ...
        'runtime_std_s' ...
    } ...
);

csvfile = fullfile(outdir, 'benchmark_topotoolbox.csv');

if isfile(csvfile)
    Told = readtable(csvfile);
    T = [Told; Tnew];
else
    T = Tnew;
end

writetable(T, csvfile);

fprintf('Saved benchmark CSV: %s\n', csvfile);

%{
% -------------------------------------------------------
% Plot DEM
% -------------------------------------------------------

fig = figure('Visible','off');

imageschs(DEM);

title('DEM');
colorbar;

saveas(fig, fullfile(outdir, 'topotoolbox_dem.png'));

close(fig);

% -------------------------------------------------------
% Plot flow accumulation
% -------------------------------------------------------

fig = figure('Visible','off');

imageschs(DEM, log10(A));

title('Flow accumulation');

colorbar;

saveas(fig, fullfile(outdir, 'topotoolbox_area.png'));

close(fig);

% -------------------------------------------------------
% Plot drainage basins
% -------------------------------------------------------

fig = figure('Visible','off');

imageschs(DEM, DB);

title('Drainage basins');

colorbar;

saveas(fig, fullfile(outdir, 'topotoolbox_catchments.png'));

close(fig);

% -------------------------------------------------------
% Hillshade + transparent basins
% -------------------------------------------------------

fig = figure('Visible','off');

imageschs(DEM);

hold on

h = imagesc(DB.Z);

set(h, 'AlphaData', 0.45);

axis image

title('Drainage basins over DEM hillshade');

colorbar;

saveas(fig, fullfile(outdir, 'topotoolbox_catchments_hillshade.png'));

close(fig);
%}
