clear; clc; close all;

%%
zipPath = 'D:\桌面\7台滑翔机数据.zip';
workDir = fullfile(tempdir, 'glider7_unpacked');

if exist(workDir, 'dir')
    rmdir(workDir, 's');
end
mkdir(workDir);
unzip(zipPath, workDir);

depthGrid  = (0:1:1000)';      
timeGridHr = -24:0.25:24;     

climRange = [-5, 5];

depthBinEdges = 0:5:1000;      
depthBinCtrs   = depthBinEdges(1:end-1) + diff(depthBinEdges)/2;

% 
sigmaZ = 1.8;   % 深度方向高斯平滑尺度（网格点）
sigmaT = 1.0;   % 时间方向高斯平滑尺度（网格点）

% 
gliderInfo = struct( ...
    'id',   {'030','MR001','039','017','002','032','010'}, ...
    't0',   { ...
        datetime(2024,11,8,23,45,0,'TimeZone','UTC'), ...
        datetime(2024,11,9, 8,39,0,'TimeZone','UTC'), ...
        datetime(2024,11,10,14,46,0,'TimeZone','UTC'), ...
        datetime(2024,11,10,15,26,0,'TimeZone','UTC'), ...
        datetime(2024,11,9, 0,23,0,'TimeZone','UTC'), ...
        datetime(2024,11,10,14,37,0,'TimeZone','UTC'), ...
        datetime(2024,11,10,18,50,0,'TimeZone','UTC') } );

% 
plotOrder = {'039','017','030','MR001','002','032','010'};

nG = numel(plotOrder);

%% 
d = dir(workDir);
d = d([d.isdir]);  
folderNames = {d.name};
folderNames = folderNames(~ismember(folderNames, {'.','..'}));

results = struct('id', [], 'folder', [], 't0', [], 'score', [], ...
                 'tGridHr', [], 'depthGrid', [], 'dcGrid', []);

for k = 1:nG
    gid = plotOrder{k};

    infoIdx = find(strcmp({gliderInfo.id}, gid), 1, 'first');
    if isempty(infoIdx)
        error('未在 gliderInfo 中找到编号 %s。', gid);
    end
    t0 = gliderInfo(infoIdx).t0;

    folderIdx = find(contains(folderNames, gid), 1, 'first');
    if isempty(folderIdx)
        error('未找到与滑翔机编号 %s 对应的文件夹。', gid);
    end
    folderName = folderNames{folderIdx};
    folderPath = fullfile(workDir, folderName);

    ncFiles = dir(fullfile(folderPath, '*.nc'));
    if isempty(ncFiles)
        error('文件夹 %s 内没有 nc 文件。', folderName);
    end


    tAll = [];
    zAll = [];
    TAll = [];
    SAll = [];
    qcTAll = [];
    qcSAll = [];
    qcZAll = [];

    for j = 1:numel(ncFiles)
        ncPath = fullfile(folderPath, ncFiles(j).name);

        rawTime = double(ncread(ncPath, 'TIME'));
        depth   = double(ncread(ncPath, 'DEPTH'));
        temp    = double(ncread(ncPath, 'TEMPERATURE'));
        sal     = double(ncread(ncPath, 'SALINITY'));

        try, qcT = double(ncread(ncPath, 'TEMPERATURE_QC')); catch, qcT = nan(size(temp)); end
        try, qcS = double(ncread(ncPath, 'SALINITY_QC'));     catch, qcS = nan(size(sal));  end
        try, qcZ = double(ncread(ncPath, 'DEPTH_QC'));        catch, qcZ = nan(size(depth)); end


        t = datetime(1970,1,1,'TimeZone','UTC') + days(rawTime);

        tAll   = [tAll; t(:)];
        zAll   = [zAll; depth(:)];
        TAll   = [TAll; temp(:)];
        SAll   = [SAll; sal(:)];
        qcTAll = [qcTAll; qcT(:)];
        qcSAll = [qcSAll; qcS(:)];
        qcZAll = [qcZAll; qcZ(:)];
    end

    % ---------------- 原始点清洗 ----------------
    valid = isfinite(zAll) & isfinite(TAll) & isfinite(SAll) & ~isnat(tAll);
    valid = valid & (zAll >= 0) & (zAll <= 1000);
    valid = valid & (TAll > -2) & (TAll < 40);
    valid = valid & (SAll > 20) & (SAll < 40);

    if any(isfinite(qcTAll)), valid = valid & (qcTAll <= 1 | isnan(qcTAll)); end
    if any(isfinite(qcSAll)), valid = valid & (qcSAll <= 1 | isnan(qcSAll)); end
    if any(isfinite(qcZAll)), valid = valid & (qcZAll <= 1 | isnan(qcZAll)); end

    tAll = tAll(valid);
    zAll = zAll(valid);
    TAll = TAll(valid);
    SAll = SAll(valid);

    % 经典声速公式
    cAll = mackenzie1981(TAll, SAll, zAll);


    cValid = isfinite(cAll) & (cAll > 1430) & (cAll < 1560);
    tAll = tAll(cValid);
    zAll = zAll(cValid);
    cAll = cAll(cValid);


    tRelHr = hours(tAll - t0);


    preMask = (tRelHr >= -24) & (tRelHr < 0);

    zPre = zAll(preMask);
    cPrePts = cAll(preMask);

    binIdx = discretize(zPre, depthBinEdges);
    cPreBin = nan(numel(depthBinCtrs), 1);

    for ib = 1:numel(depthBinCtrs)
        v = cPrePts(binIdx == ib);
        if ~isempty(v)
            cPreBin(ib) = median(v, 'omitnan');
        end
    end

    cPreBin = fillmissing(cPreBin, 'linear', 'EndValues', 'nearest');
    cPreBin = smoothdata(cPreBin, 'movmean', 5, 'omitnan');

    cPreAtRaw = interp1(depthBinCtrs, cPreBin, zAll, 'linear', 'extrap');


    dcRaw = cAll - cPreAtRaw;


    keep = isfinite(dcRaw) & ~isoutlier(dcRaw, 'median', 'ThresholdFactor', 6);
    tRelHr = tRelHr(keep);
    zAll   = zAll(keep);
    dcRaw  = dcRaw(keep);


    [TG, ZG] = meshgrid(timeGridHr, depthGrid);
    dcGrid = griddata(tRelHr, zAll, dcRaw, TG, ZG, 'linear');


    nanMask = isnan(dcGrid);
    if any(nanMask(:))
        dcGrid(nanMask) = griddata(tRelHr, zAll, dcRaw, TG(nanMask), ZG(nanMask), 'nearest');
    end


    dcGrid = gaussianSmooth2D(dcGrid, sigmaZ, sigmaT);


    dcGrid = max(min(dcGrid, 8), -8);


    upperMask = depthGrid <= 300;
    score = -mean(dcGrid(upperMask, :), 'all', 'omitnan');

    results(k).id = gid;
    results(k).folder = folderName;
    results(k).t0 = t0;
    results(k).score = score;
    results(k).tGridHr = timeGridHr;
    results(k).depthGrid = depthGrid;
    results(k).dcGrid = dcGrid;
end

%%
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 56 14]);
tl = tiledlayout(fig, 1, nG, 'TileSpacing', 'compact', 'Padding', 'compact');

cmap = blueWhiteRed(256);
colormap(fig, cmap);

for k = 1:nG
    ax = nexttile(tl, k);


    imagesc(ax, results(k).tGridHr, results(k).depthGrid, results(k).dcGrid);
    hold(ax, 'on');


    xline(ax, 0, '--', ...
        'Color', [0.85 0 0], ...
        'LineWidth', 1.5);


    try
        contour(ax, results(k).tGridHr, results(k).depthGrid, results(k).dcGrid, ...
            [0 0], 'k-', 'LineWidth', 0.8);
    catch
    end

    set(ax, ...
        'YDir', 'reverse', ...
        'LineWidth', 1.4, ...
        'FontName', 'Arial', ...
        'FontSize', 8, ...
        'TickDir', 'out', ...
        'Box', 'off', ...
        'XColor', 'k', ...
        'YColor', 'k');

    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';

    xlim(ax, [-24 24]);
    ylim(ax, [0 1000]);
    caxis(ax, climRange);


    xticks(ax, -24:12:24);
    yticks(ax, 0:200:1000);

    xticklabels(ax, {'-24','-12','0','12','24'});
    yticklabels(ax, {'0','200','400','600','800','1000'});
end


cb = colorbar;
cb.Position = [0.915 0.15 0.015 0.72];
cb.Limits = climRange;
cb.Ticks = -5:1:5;
cb.FontName = 'Arial';
cb.FontSize = 9;
cb.LineWidth = 1;
cb.Label.String = '\Delta c (m s^{-1})';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 10;

set(tl, 'Padding', 'compact', 'TileSpacing', 'compact');


outPng = fullfile(workDir, '7_gliders_sound_speed_anomaly_final_with_colorbar.png');
exportgraphics(fig, outPng, 'Resolution', 600);

disp(['图已导出到: ' outPng]);

%%

function c = mackenzie1981(T, S, D)
% Mackenzie, K. V. (1981) sound speed formula
    c = 1448.96 ...
        + 4.591 .* T ...
        - 5.304e-2 .* T.^2 ...
        + 2.374e-4 .* T.^3 ...
        + 1.340 .* (S - 35) ...
        + 1.630e-2 .* D ...
        + 1.675e-7 .* D.^2 ...
        - 1.025e-2 .* T .* (S - 35) ...
        - 7.139e-13 .* T .* D.^3;
end

function cmap = blueWhiteRed(n)

    if nargin < 1
        n = 256;
    end

    anchors = [ ...
        0.06 0.20 0.58; ...
        0.35 0.58 0.93; ...
        1.00 1.00 1.00; ...
        0.97 0.68 0.50; ...
        0.76 0.08 0.12];

    x0 = linspace(0, 1, size(anchors, 1));
    xq = linspace(0, 1, n);
    cmap = interp1(x0, anchors, xq, 'pchip');
    cmap = max(min(cmap, 1), 0);
end

function B = gaussianSmooth2D(A, sigmaZ, sigmaT)

    if nargin < 2, sigmaZ = 1.5; end
    if nargin < 3, sigmaT = 1.0; end

    halfZ = max(1, ceil(3*sigmaZ));
    halfT = max(1, ceil(3*sigmaT));

    [x, y] = meshgrid(-halfT:halfT, -halfZ:halfZ);
    G = exp(-(x.^2/(2*sigmaT^2) + y.^2/(2*sigmaZ^2)));
    G = G ./ sum(G(:));

    B = conv2(A, G, 'same');
end