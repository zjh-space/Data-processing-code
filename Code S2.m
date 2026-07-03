clear; clc; close all;

%% 
zipFile = "D:\桌面\No1-7.zip";

if ~isfile(zipFile)
    cand = dir("D:\桌面\No1-7*.zip");
    if isempty(cand)
        error("找不到 zip 文件，请检查路径。");
    end
    zipFile = fullfile(cand(1).folder, cand(1).name);
end

extractDir = fullfile(tempdir, "glider_fig4_extract");
if exist(extractDir, "dir")
    rmdir(extractDir, "s");
end
mkdir(extractDir);

outPng = fullfile(fileparts(zipFile), "Fig4_sound_speed_structure.png");
outPdf = fullfile(fileparts(zipFile), "Fig4_sound_speed_structure.pdf");

% Common depth grid
zGrid = (0:1:780).';   % m

% Smoothing for display only
smoothWindow = 7;

% Inset zoom depth range
zoomDepthMax = 220;

%% 
unzip(zipFile, extractDir);

%% 
d = dir(fullfile(extractDir, "**", "*.nc"));
if isempty(d)
    error("未找到 nc 文件。");
end
fprintf("找到 %d 个 nc 文件\n", numel(d));

%%
topFolders = strings(numel(d), 1);
for i = 1:numel(d)
    relPath = string(d(i).folder);
    relPath = erase(relPath, string(extractDir) + filesep);
    parts = split(relPath, filesep);
    if isempty(parts) || strlength(parts(1)) == 0
        topFolders(i) = "root";
    else
        topFolders(i) = parts(1);
    end
end

gliderFolders = unique(topFolders, "stable");
gliderFolders = sort(gliderFolders);
nGlider = numel(gliderFolders);

fprintf("识别出 %d 台滑翔机文件夹\n", nGlider);

%%
gliderLabel = cell(nGlider, 1);

gliderMeanC    = NaN(numel(zGrid), nGlider);
gliderStdC     = NaN(numel(zGrid), nGlider);
gliderMeanGrad = NaN(numel(zGrid), nGlider);
gliderStdGrad  = NaN(numel(zGrid), nGlider);

upperMeanC   = NaN(nGlider, 1);
upperRangeC  = NaN(nGlider, 1);
maxGradDepth = NaN(nGlider, 1);
minCDepth    = NaN(nGlider, 1);
profileCount = zeros(nGlider, 1);

%% 
for g = 1:nGlider
    folderName = gliderFolders(g);
    gliderLabel{g} = normalizeLabel(folderName);

    gliderPath = fullfile(extractDir, folderName);
    f = dir(fullfile(gliderPath, "**", "*.nc"));

    if isempty(f)
        warning("文件夹中没有 nc 文件: %s", gliderPath);
        continue;
    end

    cStack = [];
    gradStack = [];

    for k = 1:numel(f)
        filePath = fullfile(f(k).folder, f(k).name);

        try
            profiles = readProfilesFromNc(filePath);

            if isempty(profiles)
                continue;
            end

            for p = 1:numel(profiles)
                z = profiles{p}.z(:);
                T = profiles{p}.T(:);
                S = profiles{p}.S(:);

                valid = isfinite(z) & isfinite(T) & isfinite(S);
                z = z(valid);
                T = T(valid);
                S = S(valid);

                if numel(z) < 20
                    continue;
                end

                [z, idxSort] = sort(z);
                T = T(idxSort);
                S = S(idxSort);

                [z, idxUniq] = unique(z, "stable");
                T = T(idxUniq);
                S = S(idxUniq);

                if numel(z) < 20
                    continue;
                end

                c = soundSpeedMackenzie(T, S, z);
                cI = interp1(z, c, zGrid, "pchip", NaN);

                if sum(isfinite(cI)) < 50
                    continue;
                end

                cStack(:, end+1) = cI; %#ok<SAGROW>
                gradStack(:, end+1) = gradient(cI, zGrid); %#ok<SAGROW>
                profileCount(g) = profileCount(g) + 1;
            end

        catch ME
            warning("%s 读取失败: %s", filePath, ME.message);
        end
    end

    if isempty(cStack)
        warning("该文件夹没有保留到有效剖面: %s", folderName);
        continue;
    end

    gliderMeanC(:, g)    = mean(cStack, 2, "omitnan");
    gliderStdC(:, g)     = std(cStack, 0, 2, "omitnan");
    gliderMeanGrad(:, g) = mean(gradStack, 2, "omitnan");
    gliderStdGrad(:, g)  = std(gradStack, 0, 2, "omitnan");

    gliderMeanC(:, g)    = smoothdata(gliderMeanC(:, g), "movmean", smoothWindow);
    gliderMeanGrad(:, g) = smoothdata(gliderMeanGrad(:, g), "movmean", smoothWindow);

    cThis = gliderMeanC(:, g);
    gradThis = gliderMeanGrad(:, g);

    idx100 = zGrid <= 100;
    idx200 = zGrid <= 200;

    upperMeanC(g) = mean(cThis(idx100), "omitnan");

    c200 = cThis(idx200);
    upperRangeC(g) = max(c200, [], "omitnan") - min(c200, [], "omitnan");

    grad200 = gradThis(idx200);
    [~, ii] = max(abs(grad200));
    z200 = zGrid(idx200);
    maxGradDepth(g) = z200(ii);

    [~, jj] = min(cThis);
    minCDepth(g) = zGrid(jj);
end

%%
validGlider = any(isfinite(gliderMeanC), 1);

gliderLabel    = gliderLabel(validGlider);
gliderMeanC    = gliderMeanC(:, validGlider);
gliderStdC     = gliderStdC(:, validGlider);
gliderMeanGrad = gliderMeanGrad(:, validGlider);
gliderStdGrad  = gliderStdGrad(:, validGlider);

upperMeanC   = upperMeanC(validGlider);
upperRangeC  = upperRangeC(validGlider);
maxGradDepth = maxGradDepth(validGlider);
minCDepth    = minCDepth(validGlider);
profileCount = profileCount(validGlider);

nGlider = numel(gliderLabel);

if nGlider < 2
    error("有效滑翔机数量过少，无法绘图。");
end

fprintf("最终保留 %d 台滑翔机\n", nGlider);

%% 
ensembleMeanC    = mean(gliderMeanC, 2, "omitnan");
ensembleStdC     = std(gliderMeanC, 0, 2, "omitnan");
ensembleMeanGrad = mean(gliderMeanGrad, 2, "omitnan");
ensembleStdGrad  = std(gliderMeanGrad, 0, 2, "omitnan");

ensembleMeanC    = smoothdata(ensembleMeanC, "movmean", smoothWindow);
ensembleStdC     = smoothdata(ensembleStdC, "movmean", smoothWindow);
ensembleMeanGrad = smoothdata(ensembleMeanGrad, "movmean", smoothWindow);
ensembleStdGrad  = smoothdata(ensembleStdGrad, "movmean", smoothWindow);

gliderAnomC = gliderMeanC - ensembleMeanC;

%% 
fig = figure("Color", "w", "Position", [80, 60, 1700, 1200]);
set(fig, "Renderer", "painters");

t = tiledlayout(2, 2, "Padding", "compact", "TileSpacing", "compact");

colors = lines(nGlider);

applyAxisStyle = @(ax) set(ax, ...
    "FontName", "Times New Roman", ...
    "FontSize", 11, ...
    "LineWidth", 1.6, ...
    "Box", "off", ...
    "TickDir", "out", ...
    "TickLength", [0.012 0.012], ...
    "XMinorTick", "off", ...
    "YMinorTick", "off");

axesForFrame = gobjects(0);

%% 
ax1 = nexttile(t, 1);
hold(ax1, "on");

fillBand(ax1, ensembleMeanC, ensembleStdC, zGrid, [0.2 0.2 0.2], 0.12);

hGl = gobjects(nGlider, 1);
for g = 1:nGlider
    hGl(g) = plot(ax1, gliderMeanC(:, g), zGrid, ...
        "Color", [colors(g, :) 0.90], ...
        "LineWidth", 1.25);
end
hEns = plot(ax1, ensembleMeanC, zGrid, "k-", "LineWidth", 2.8);

set(ax1, "YDir", "reverse");
ylim(ax1, [0 max(zGrid)]);
xlabel(ax1, "Sound speed / m s^{-1}");
ylabel(ax1, "Depth / m");
title(ax1, "Mean sound-speed profiles", "FontWeight", "bold");
xlim(ax1, [min(gliderMeanC, [], "all", "omitnan") - 1.2, max(gliderMeanC, [], "all", "omitnan") + 1.2]);
applyAxisStyle(ax1);

lgd = legend(ax1, [hEns; hGl], ["Ensemble mean"; string(gliderLabel(:))], ...
    "Location", "northeast", "FontSize", 7.5, "Box", "off");
lgd.NumColumns = 2;

% 
drawnow;
mainPos = ax1.Position;
insetPos = [mainPos(1) + 0.48 * mainPos(3), ...
            mainPos(2) + 0.08 * mainPos(4), ...
            0.33 * mainPos(3), ...
            0.24 * mainPos(4)];

ax1Inset = axes("Position", insetPos, ...
    "Box", "off", ...
    "Color", "w", ...
    "FontName", "Times New Roman", ...
    "FontSize", 8, ...
    "LineWidth", 1.2, ...
    "TickDir", "out", ...
    "TickLength", [0.010 0.010], ...
    "XMinorTick", "off", ...
    "YMinorTick", "off");
hold(ax1Inset, "on");

idxZoom = zGrid <= zoomDepthMax;
fillBand(ax1Inset, ensembleMeanC(idxZoom), ensembleStdC(idxZoom), zGrid(idxZoom), [0.2 0.2 0.2], 0.18);
plot(ax1Inset, ensembleMeanC(idxZoom), zGrid(idxZoom), "k-", "LineWidth", 2.4);

set(ax1Inset, "YDir", "reverse");
ylim(ax1Inset, [0 zoomDepthMax]);

xZoomMin = min(ensembleMeanC(idxZoom) - ensembleStdC(idxZoom), [], "omitnan");
xZoomMax = max(ensembleMeanC(idxZoom) + ensembleStdC(idxZoom), [], "omitnan");
xMargin = max(0.8, 0.02 * (xZoomMax - xZoomMin));
xlim(ax1Inset, [xZoomMin - xMargin, xZoomMax + xMargin]);

insetXLim = get(ax1Inset, "XLim");
xticks(ax1Inset, round(linspace(insetXLim(1), insetXLim(2), 4), 0));
yticks(ax1Inset, [0 50 100 150 200]);
xlabel(ax1Inset, "c", "FontSize", 8);
ylabel(ax1Inset, "z / m", "FontSize", 8);
applyAxisStyle(ax1Inset);

axesForFrame(end+1) = ax1; 
axesForFrame(end+1) = ax1Inset; 

%% 
ax2 = nexttile(t, 2);
hold(ax2, "on");

fillBand(ax2, zeros(size(ensembleMeanC)), ensembleStdC, zGrid, [0.2 0.2 0.2], 0.10);
for g = 1:nGlider
    plot(ax2, gliderAnomC(:, g), zGrid, ...
        "Color", [colors(g, :) 0.90], ...
        "LineWidth", 1.25);
end
xline(ax2, 0, "k--", "LineWidth", 1.0);

set(ax2, "YDir", "reverse");
ylim(ax2, [0 max(zGrid)]);
xlabel(ax2, "\Delta c / m s^{-1}");
ylabel(ax2, "Depth / m");
title(ax2, "Sound-speed anomaly relative to ensemble mean", "FontWeight", "bold");
xlim(ax2, [min(gliderAnomC, [], "all", "omitnan") - 0.3, max(gliderAnomC, [], "all", "omitnan") + 0.3]);
applyAxisStyle(ax2);

axesForFrame(end+1) = ax2; 

%%
ax3 = nexttile(t, 3);
hold(ax3, "on");

fillBand(ax3, ensembleMeanGrad * 1e3, ensembleStdGrad * 1e3, zGrid, [0.2 0.2 0.2], 0.12);
for g = 1:nGlider
    plot(ax3, gliderMeanGrad(:, g) * 1e3, zGrid, ...
        "Color", [colors(g, :) 0.90], ...
        "LineWidth", 1.25);
end
xline(ax3, 0, "k--", "LineWidth", 1.0);

set(ax3, "YDir", "reverse");
ylim(ax3, [0 max(zGrid)]);
xlabel(ax3, "Vertical sound-speed gradient (10^{-3} s^{-1})");
ylabel(ax3, "Depth / m");
title(ax3, "Sound-speed vertical gradient", "FontWeight", "bold");
xlim(ax3, [min(gliderMeanGrad, [], "all", "omitnan") * 1e3 - 1, max(gliderMeanGrad, [], "all", "omitnan") * 1e3 + 1]);
applyAxisStyle(ax3);

axesForFrame(end+1) = ax3; 

%% 
ax4 = nexttile(t, 4);
hold(ax4, "on");
layerEdges = [0 50 150 300 780];
layerNames = { ...
    "0–50 m", ...
    "50–150 m", ...
    "150–300 m", ...
    "300–780 m"};

nLayer = numel(layerEdges) - 1;
layerContribution = NaN(nGlider, nLayer);

for g = 1:nGlider
    anom = abs(gliderAnomC(:, g));
    anom(~isfinite(anom)) = 0;

    totalA = trapz(zGrid, anom);
    if totalA <= 0
        continue;
    end

    raw = zeros(1, nLayer);
    for L = 1:nLayer
        if L < nLayer
            idx = zGrid >= layerEdges(L) & zGrid < layerEdges(L+1);
        else
            idx = zGrid >= layerEdges(L) & zGrid <= layerEdges(L+1);
        end

        if any(idx)
            raw(L) = trapz(zGrid(idx), anom(idx));
        else
            raw(L) = 0;
        end
    end

    if sum(raw) > 0
        row = 100 * raw / sum(raw);
        row = row / sum(row) * 100;    
        row(4) = 100 - sum(row(1:3)); 
        layerContribution(g, :) = row;
    end
end

% 
desiredOrder = ["030","001","039","017","002","032","010"];
currentLabels = string(gliderLabel(:));

orderIdx = [];
for ii = 1:numel(desiredOrder)
    k = find(currentLabels == desiredOrder(ii), 1);
    if ~isempty(k)
        orderIdx(end+1) = k; 
    else
        warning("Fig.4d order label not found: %s", desiredOrder(ii));
    end
end


remainingIdx = setdiff(1:nGlider, orderIdx, "stable");
orderIdx = [orderIdx(:); remainingIdx(:)];

layerContributionOrd = layerContribution(orderIdx, :);
labelOrd = gliderLabel(orderIdx);


b = bar(ax4, layerContributionOrd, ...
    "stacked", ...
    "BarWidth", 0.75, ...
    "LineWidth", 1.0, ...
    "EdgeColor", "none");


layerColors = [
    0.98 0.78 0.95
    0.94 0.56 0.89
    0.86 0.30 0.78
    0.68 0.10 0.62
];

for L = 1:nLayer
    b(L).FaceColor = layerColors(L, :);
end

xticks(ax4, 1:nGlider);
xticklabels(ax4, labelOrd);
xtickangle(ax4, 0);

ylim(ax4, [0 100]);
xlim(ax4, [0.4, nGlider + 0.6]);
xlabel(ax4, "Glider");
ylabel(ax4, "Contribution / %");
title(ax4, "Depth-layer contribution to sound-speed anomaly", "FontWeight", "bold");

legend(ax4, layerNames, ...
    "Location", "northeast", ...
    "Box", "off", ...
    "FontSize", 8);

applyAxisStyle(ax4);


for j = 1:nGlider
    if all(isfinite(layerContributionOrd(j, :)))
        totalTxt = sprintf("100%%");
        text(ax4, j, 101.2, totalTxt, ...
            "HorizontalAlignment", "center", ...
            "VerticalAlignment", "bottom", ...
            "FontSize", 7.5, ...
            "FontWeight", "bold");
    end
end

axesForFrame(end+1) = ax4; 

sgtitle(t, "Fig.4  Sound-speed structure reconstruction from seven gliders", ...
    "FontName", "Times New Roman", "FontWeight", "bold", "FontSize", 15);

%% 
drawnow;
for k = 1:numel(axesForFrame)
    drawTopRightFrame(fig, axesForFrame(k), 1.6);
end

%%
% exportgraphics(fig, outPng, "Resolution", 600);
% exportgraphics(fig, outPdf, "ContentType", "vector");
% 
% fprintf("Fig.4 已保存到:\n%s\n%s\n", outPng, outPdf);

%%
function profiles = readProfilesFromNc(filePath)

    info = ncinfo(filePath);
    varNames = string({info.Variables.Name});

    Tname = pickVarName(varNames, ["TEMPERATURE","TEMP","T"]);
    Sname = pickVarName(varNames, ["SALINITY","PSAL","SAL"]);
    Zname = pickVarName(varNames, ["DEPTH","PRESSURE","PRES","P","Z"]);

    Traw = double(ncread(filePath, char(Tname)));
    Sraw = double(ncread(filePath, char(Sname)));
    Zraw = double(ncread(filePath, char(Zname)));

    Traw = squeeze(Traw);
    Sraw = squeeze(Sraw);
    Zraw = squeeze(Zraw);

    profiles = {};

    if isvector(Traw) && isvector(Sraw) && isvector(Zraw)
        profiles{1} = struct("z", Zraw(:), "T", Traw(:), "S", Sraw(:));
        return;
    end

    [Z, T, S] = orientToProfiles(Zraw, Traw, Sraw);

    if isvector(T) && isvector(S) && isvector(Z)
        profiles{1} = struct("z", Z(:), "T", T(:), "S", S(:));
        return;
    end

    if ismatrix(T) && ismatrix(S) && ismatrix(Z)
        nProf = min([size(T,2), size(S,2), size(Z,2)]);
        profiles = cell(nProf, 1);
        for p = 1:nProf
            profiles{p} = struct("z", Z(:, p), "T", T(:, p), "S", S(:, p));
        end
        return;
    end

    profiles{1} = struct("z", Zraw(:), "T", Traw(:), "S", Sraw(:));
end

function [Z, T, S] = orientToProfiles(Z, T, S)
% 

    Z = squeeze(Z);
    T = squeeze(T);
    S = squeeze(S);

    if isvector(Z), Z = Z(:); end
    if isvector(T), T = T(:); end
    if isvector(S), S = S(:); end

    if ismatrix(T) && ismatrix(S) && ismatrix(Z)
        if isequal(size(T), size(S), size(Z))
            return;
        end

        if isequal(size(T.'), size(S.')) && isequal(size(T.'), size(Z.'))
            Z = Z.'; T = T.'; S = S.';
            return;
        end

        if size(T,1) == size(S,1)
            if size(Z,1) ~= size(T,1) && size(Z,2) == size(T,1)
                Z = Z.';
            end
            return;
        end

        if size(T,2) == size(S,2)
            if size(Z,2) ~= size(T,2) && size(Z,1) == size(T,2)
                Z = Z.';
            end
            T = T.'; S = S.';
            return;
        end
    end

    Z = Z(:);
    T = T(:);
    S = S(:);
end

function vName = pickVarName(varNames, candidates)


    vName = "";
    lowerNames = lower(varNames);

    for c = candidates
        idx = find(contains(lowerNames, lower(c)), 1);
        if ~isempty(idx)
            vName = varNames(idx);
            return;
        end
    end

    error("变量未找到，尝试的关键词是: %s", strjoin(string(candidates), ", "));
end

function c = soundSpeedMackenzie(T, S, z)


    T = double(T);
    S = double(S);
    z = double(z);

    c = 1448.96 ...
        + 4.591*T ...
        - 5.304e-2*T.^2 ...
        + 2.374e-4*T.^3 ...
        + 1.340*(S - 35) ...
        + 1.630e-2*z ...
        + 1.675e-7*z.^2 ...
        - 1.025e-2*T.*(S - 35) ...
        - 7.139e-13*T.*z.^3;
end

function fillBand(ax, centerLine, spreadLine, z, rgb, faceAlpha)


    centerLine = centerLine(:);
    spreadLine = spreadLine(:);
    z = z(:);

    valid = isfinite(centerLine) & isfinite(spreadLine) & isfinite(z);
    centerLine = centerLine(valid);
    spreadLine = spreadLine(valid);
    z = z(valid);

    if isempty(centerLine)
        return;
    end

    x1 = centerLine - spreadLine;
    x2 = centerLine + spreadLine;

    xx = [x1; flipud(x2)];
    yy = [z; flipud(z)];

    patch(ax, xx, yy, rgb, ...
        "FaceAlpha", faceAlpha, ...
        "EdgeColor", "none");
end

function cmap = redBlueCmap(n)

    if nargin < 1
        n = 256;
    end
    n1 = floor(n/2);
    n2 = n - n1;

    blue  = [0.10 0.25 0.80];
    white = [1 1 1];
    red   = [0.85 0.15 0.10];

    c1 = [linspace(blue(1), white(1), n1)', ...
          linspace(blue(2), white(2), n1)', ...
          linspace(blue(3), white(3), n1)'];

    c2 = [linspace(white(1), red(1), n2)', ...
          linspace(white(2), red(2), n2)', ...
          linspace(white(3), red(3), n2)'];

    cmap = [c1; c2];
end

function label = normalizeLabel(folderName)

    folderName = string(folderName);
    digits = regexp(folderName, "\d+", "match", "once");
    if ~isempty(digits)
        label = char(digits);
    else
        label = char(folderName);
    end
end

function drawTopRightFrame(fig, ax, lw)

    if ~isvalid(ax)
        return;
    end

    pos = ax.Position;
    x = pos(1);
    y = pos(2);
    w = pos(3);
    h = pos(4);

    annotation(fig, "line", [x, x+w], [y+h, y+h], ...
        "Color", "k", "LineWidth", lw);
    annotation(fig, "line", [x+w, x+w], [y, y+h], ...
        "Color", "k", "LineWidth", lw);
end