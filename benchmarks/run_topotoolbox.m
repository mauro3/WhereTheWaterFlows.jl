clear; clc;

%We compare the time required to compute single-flow-direction routing, flow accumulation, and drainage basin labels from the same DEM. 
% For WWF this is one waterflows call; for TopoToolbox this is FLOWobj + flowacc + drainagebasins.

maxNumCompThreads(1);
nthreads = maxNumCompThreads;

fprintf('MATLAB threads: %d\n', nthreads);

% Reproducible paths
scriptdir = fileparts(mfilename('fullpath'));
root = fileparts(scriptdir);

topodir = fullfile(fileparts(root), 'topotoolbox3');

if ~exist(topodir, 'dir')
    error('TopoToolbox not found. Clone it next to this repository: git clone https://github.com/TopoToolbox/topotoolbox3.git')
end

addpath(genpath(topodir))

datadir = fullfile(root, 'data_raw');
outdir  = fullfile(root, 'outputs', 'topotoolbox');

demfile_small = fullfile(datadir, 'swissalti3d_tile.tif');
demfile_large = fullfile(datadir, 'swissalti3d_tilelarge.tif');

if ~exist(demfile_large, 'file')
    fprintf('Creating large DEM for TopoToolbox...\n');

    [Z, R] = readgeoraster(demfile_small);
    Z = single(Z);

    Zlarge = imresize(Z, 4, 'bicubic');

    Rlarge = R;
    Rlarge.RasterSize = size(Zlarge);
    Rlarge.CellExtentInWorldX = R.CellExtentInWorldX / 4;
    Rlarge.CellExtentInWorldY = R.CellExtentInWorldY / 4;

    geotiffwrite(demfile_large, Zlarge, Rlarge, ...
        'CoordRefSysCode', 2056);

    fprintf('Saved large DEM to: %s\n', demfile_large);
end

demfile = demfile_large; % or change to demfile_small


DEM = GRIDobj(demfile);

% -------------------------------------------------------
% Load DEM
% -------------------------------------------------------

DEM = GRIDobj(demfile);

npixels = nnz(~isnan(DEM.Z));


% -------------------------------------------------------
% Run TopoToolbox benchmark
% -------------------------------------------------------

nruns = 6;

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

mean_runtime = mean(runtime_total_s(3:end));
std_runtime  = std(runtime_total_s(3:end));

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



