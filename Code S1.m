%% Typhoon-induced acoustic propagation from glider CTD data
% MATLAB R2024b
%
% Note:
%   The TL solver below is a compact MMPE-style split-step parabolic-equation
%   implementation intended for figure production and sensitivity comparison.
%   For formal operational prediction, replace runMMPEstylePE() with a
%   validated MMPE/RAM/BELLHOP executable using the same SSPs.

clear; close all; clc;


zipFile = "D:\桌面\No8.zip";
outDir  = fullfile(pwd, "glider_acoustic_output");

typhoonStart = datetime(2024,11,8,9,0,0);
typhoonEnd   = datetime(2024,11,11,10,0,0);

sourceDepthM   = 50;     % typical acoustic source depth
receiverDepthM = 50;     % panel (c) receiver depth
freqHz         = 400;    % acoustic frequency for PE/TL calculation
maxRangeKm     = 50;     % horizontal range shown in panels

zBinM       = 2;         % CTD bin size for SSP construction
peDzM       = 2;         % PE depth step
peDrM       = 50;        % PE range step
bottomPadM  = 40;        % extra depth below deepest CTD bin for sponge layer
rayAnglesDeg = [-4 -2 0 2 4]; % a few representative launch angles, kept sparse for clarity
nBootSSP = 120;       % bootstrap members for SSP uncertainty
nBootTL  = 18;        % subset propagated through PE for TL uncertainty
rng(42);              % reproducible uncertainty envelopes

if ~exist(outDir, "dir")
    mkdir(outDir);
end

%% 
dataRoot = fullfile(outDir, "unzipped_1hao");
if ~exist(dataRoot, "dir") || isempty(dir(fullfile(dataRoot, "**", "*.nc")))
    fprintf("Unzipping %s ...\n", zipFile);
    unzip(zipFile, dataRoot);
end

files = dir(fullfile(dataRoot, "**", "*_ctd.nc"));
assert(~isempty(files), "No *_ctd.nc files found after unzipping %s.", zipFile);

G = readGliderCTD(files);
G = sortrows(G, "Time");

maxObsDepth = ceil(max(G.Depth, [], "omitnan")/10)*10;
zProfile = (0:zBinM:maxObsDepth).';

periods = struct( ...
    "name",  {"Before", "During", "After"}, ...
    "label", {"Pre-typhoon", "Typhoon passage", "Post-typhoon"}, ...
    "t0",    {min(G.Time), typhoonStart, typhoonEnd}, ...
    "t1",    {typhoonStart, typhoonEnd, max(G.Time)}, ...
    "color", {[0.09 0.28 0.55], [0.73 0.20 0.16], [0.10 0.45 0.33]});

profiles = repmat(makeSSP(G, periods(1), zProfile), 1, numel(periods));
for k = 2:numel(periods)
    profiles(k) = makeSSP(G, periods(k), zProfile);
end

fprintf("Bootstrapping SSP uncertainty by NetCDF profile segment (%d members) ...\n", nBootSSP);
for k = 1:numel(periods)
    cBoot = bootstrapSSP(G, periods(k), zProfile, nBootSSP);
    profiles(k).cBoot = cBoot;
    profiles(k).cLo = prctile(cBoot, 5, 2);
    profiles(k).cHi = prctile(cBoot, 95, 2);
end

fprintf("\nConstructed sound-speed profiles:\n");
for k = 1:numel(profiles)
    fprintf("  %-7s: %s to %s, n = %d, c(0 m) = %.2f m/s, c(%d m) = %.2f m/s\n", ...
        periods(k).name, string(periods(k).t0), string(periods(k).t1), ...
        profiles(k).nObs, profiles(k).c(1), round(sourceDepthM), ...
        interp1(profiles(k).z, profiles(k).c, sourceDepthM, "linear", "extrap"));
end

%% 
rangeM = 0:peDrM:(maxRangeKm*1000);
bottomM = maxObsDepth + bottomPadM;

for k = 1:numel(profiles)
    profiles(k).rays = traceRaysSnell(profiles(k).z, profiles(k).c, ...
        sourceDepthM, rayAnglesDeg, maxRangeKm*1000, bottomM);

    [profiles(k).tl, profiles(k).zPE, profiles(k).rPE] = runMMPEstylePE( ...
        profiles(k).z, profiles(k).c, sourceDepthM, freqHz, rangeM, peDzM, bottomM);

    profiles(k).tlAtReceiver = interp1(profiles(k).zPE, profiles(k).tl, ...
        receiverDepthM, "linear", "extrap");
end

fprintf("Propagating bootstrap SSP subset through PE (%d members per period) ...\n", nBootTL);
for k = 1:numel(profiles)
    bootPick = round(linspace(1, nBootSSP, nBootTL));
    tlBoot = nan(nBootTL, numel(rangeM));
    for b = 1:nBootTL
        [tlB, zB, ~] = runMMPEstylePE(profiles(k).z, profiles(k).cBoot(:,bootPick(b)), ...
            sourceDepthM, freqHz, rangeM, peDzM, bottomM);
        tlBoot(b,:) = interp1(zB, tlB, receiverDepthM, "linear", "extrap");
    end
    profiles(k).tlBootAtReceiver = smoothdata(tlBoot, 2, "movmean", 9);
    profiles(k).tlLo = prctile(profiles(k).tlBootAtReceiver, 5, 1);
    profiles(k).tlHi = prctile(profiles(k).tlBootAtReceiver, 95, 1);
end

%%
set(groot, ...
    "defaultFigureColor", "w", ...
    "defaultAxesFontName", "Arial", ...
    "defaultTextFontName", "Arial", ...
    "defaultAxesFontSize", 20, ...     
    "defaultAxesLineWidth", 2, ...   
    "defaultLineLineWidth", 1.8, ...
    "defaultColorbarFontSize", 15, ...  
    "defaultLegendFontSize", 15);       

fig = figure("Units", "centimeters", "Position", [2 1 18.3 24.5], "Color", "w");
tlayout = tiledlayout(fig, 3, 3, "TileSpacing", "compact", "Padding", "compact");

sspCLim = [min(cat(1, profiles.c), [], "omitnan"), max(cat(1, profiles.c), [], "omitnan")];
tlCLim = [45 105];
tlLevels = [55 70 85 100];
mapSSP = natureBlueRed(256);
mapTL = natureThermal(256);
rayColors = [ ...
    0.08 0.22 0.48
    0.16 0.43 0.65
    0.12 0.12 0.12
    0.73 0.35 0.18
    0.58 0.10 0.10];

% 
rayAx = gobjects(1,3);
for k = 1:numel(profiles)
    ax = nexttile(tlayout, 0 + k);
    rayAx(k) = ax;
    hold(ax, "on");
    ax.LineWidth = 1.2;  

    rPlotKm = linspace(0, maxRangeKm, 220);
    Cbg = repmat(profiles(k).c(:), 1, numel(rPlotKm));
    imagesc(ax, rPlotKm, profiles(k).z, Cbg);
    set(ax, "YDir", "reverse");
    colormap(ax, mapSSP);
    clim(ax, sspCLim);

    for j = 1:numel(profiles(k).rays)
        ray = profiles(k).rays(j);
        plot(ax, ray.r/1000, ray.z, "Color", rayColors(j,:), "LineWidth", 1.45);
    end
    plot(ax, 0, sourceDepthM, "o", "MarkerFaceColor", [0 0 0], ...
        "MarkerEdgeColor", "w", "MarkerSize", 4.8);

    %
    upper = profiles(k).z <= 250;
    [~, izGrad] = max(abs(gradient(profiles(k).c(upper), profiles(k).z(upper))));
    zUpper = profiles(k).z(upper);
    zGrad = zUpper(izGrad);
    yline(ax, zGrad, ":", "Color", [0.05 0.05 0.05], "LineWidth", 1.0);

    title(ax, sprintf("%s\n%s", periods(k).label, dateSpanText(periods(k).t0, periods(k).t1)), ...
        "FontWeight", "normal", "FontSize", 12);
    xlabel(ax, "Range (km)", "FontSize", 11);
    ylabel(ax, "Depth (m)", "FontSize", 11);
    ylim(ax, [0 bottomM]);
    xlim(ax, [0 maxRangeKm]);
    box(ax, "on");
    grid(ax, "off");
    ax.Layer = "top";
end
cb1 = colorbar(rayAx(3), "eastoutside");
cb1.Label.String = "Sound speed (m s^{-1})";
cb1.Label.FontSize = 11;
cb1.LineWidth = 1.2;

% 
tlAx = gobjects(1,3);
for k = 1:numel(profiles)
    ax = nexttile(tlayout, 3 + k);
    tlAx(k) = ax;
    hold(ax, "on");
    ax.LineWidth = 1.2; 

    tlPlot = smoothdata(profiles(k).tl, 1, "movmean", 5);
    tlPlot = smoothdata(tlPlot, 2, "movmean", 13);

    imagesc(ax, profiles(k).rPE/1000, profiles(k).zPE, tlPlot);
    set(ax, "YDir", "reverse");
    colormap(ax, mapTL);
    clim(ax, tlCLim);
    [cWhite, hWhite] = contour(ax, profiles(k).rPE/1000, profiles(k).zPE, tlPlot, ...
        tlLevels, "Color", [0.97 0.97 0.97], "LineWidth", 0.55);
    clabel(cWhite, hWhite, "FontSize", 9, "Color", [0.97 0.97 0.97], ...
        "LabelSpacing", 550);
    [cBlack, hBlack] = contour(ax, profiles(k).rPE/1000, profiles(k).zPE, tlPlot, ...
        [70 70], "Color", [0.05 0.05 0.05], "LineWidth", 0.85);
    clabel(cBlack, hBlack, "FontSize", 9, "Color", [0.05 0.05 0.05], ...
        "LabelSpacing", 700);
    plot(ax, [0 maxRangeKm], [receiverDepthM receiverDepthM], "--", ...
        "Color", [1 1 1]*0.08, "LineWidth", 0.8);

    title(ax, periods(k).label, "FontWeight", "normal", "FontSize", 12);
    xlabel(ax, "Range (km)", "FontSize", 11);
    ylabel(ax, "Depth (m)", "FontSize", 11);
    ylim(ax, [0 bottomM]);
    xlim(ax, [0 maxRangeKm]);
    box(ax, "on");
    grid(ax, "off");
    ax.Layer = "top";
end
cb2 = colorbar(tlAx(3), "eastoutside");
cb2.Label.String = "Transmission loss (dB)";
cb2.Label.FontSize = 11;
cb2.LineWidth = 1.2;

% 
ax = nexttile(tlayout, 7, [1 3]);
hold(ax, "on");
ax.LineWidth = 1.2; 

lineH = gobjects(1,3);
for k = 1:numel(profiles)
    fill(ax, [profiles(k).rPE/1000, fliplr(profiles(k).rPE/1000)], ...
        [profiles(k).tlLo, fliplr(profiles(k).tlHi)], periods(k).color, ...
        "FaceAlpha", 0.13, "EdgeColor", "none");
end
for k = 1:numel(profiles)
    lineH(k) = plot(ax, profiles(k).rPE/1000, profiles(k).tlAtReceiver, ...
        "Color", periods(k).color, "LineWidth", 1.8);
end

ylabel(ax, sprintf("TL at %g m (dB)", receiverDepthM), "FontSize", 11);
xlabel(ax, "Range (km)", "FontSize", 11);
xlim(ax, [0 maxRangeKm]);
ylim(ax, paddedLimits([cat(2, profiles.tlAtReceiver), cat(2, profiles.tlLo), cat(2, profiles.tlHi)], 4));
ax.YColor = [0.15 0.15 0.15];
grid(ax, "on");
ax.GridAlpha = 0.14;
box(ax, "on");

legend(ax, lineH, ["Pre-typhoon", "Typhoon passage", "Post-typhoon"], ...
    "Location", "northwest", "NumColumns", 3, "Box", "off", "FontSize", 10);

title(tlayout, "Glider-constrained acoustic response to typhoon passage", ...
    "FontName", "Arial", "FontSize", 13, "FontWeight", "bold");

%% 
pngFile = fullfile(outDir, "glider_typhoon_acoustic_panels.png");
pdfFile = fullfile(outDir, "glider_typhoon_acoustic_panels.pdf");
noteFile = fullfile(outDir, "glider_typhoon_acoustic_parameters_uncertainty.txt");
% exportgraphics(fig, pngFile, "Resolution", 600);
% exportgraphics(fig, pdfFile, "ContentType", "vector");
writelines([
    "Glider-constrained acoustic response to typhoon passage"
    "Data: D:\桌面\No8.zip"
    sprintf("Typhoon-passage window: %s to %s", string(typhoonStart), string(typhoonEnd))
    sprintf("Source depth: %.0f m", sourceDepthM)
    sprintf("Receiver depth: %.0f m", receiverDepthM)
    sprintf("Frequency: %.0f Hz", freqHz)
    sprintf("Maximum range: %.0f km", maxRangeKm)
    sprintf("PE range step: %.0f m", peDrM)
    sprintf("PE depth step: %.0f m", peDzM)
    sprintf("Model bottom: %.0f m with lower sponge layer", bottomM)
    sprintf("SSP bootstrap members: %d", nBootSSP)
    sprintf("TL bootstrap PE members per period: %d", nBootTL)
    "Bootstrap unit: NetCDF profile segment/file within each period"
    "Uncertainty included: period-internal CTD sampling variability propagated to SSP and receiver-depth TL"
    "Uncertainty not included: seabed geoacoustics, bathymetry, sea-surface roughness, 3-D variability, acoustic source/receiver calibration, and CTD systematic bias"
    ], noteFile);

fprintf("\nSaved figure and notes:\n  %s\n  %s\n  %s\n", pngFile, pdfFile, noteFile);

%% 
function G = readGliderCTD(files)
    Time = datetime.empty(0,1);
    Depth = [];
    Temp = [];
    Salt = [];
    Lat = [];
    Lon = [];
    FileID = [];

    for i = 1:numel(files)
        f = fullfile(files(i).folder, files(i).name);
        t = double(ncread(f, "TIME"));
        z = double(ncread(f, "DEPTH"));
        T = double(ncread(f, "TEMPERATURE"));
        S = double(ncread(f, "SALINITY"));
        la = double(ncread(f, "LATITUDE"));
        lo = double(ncread(f, "LONGITUDE"));

        tdt = decodeGliderTime(t);
        z(z < -1000) = NaN;
        T(T < -1000 | T < -3 | T > 45) = NaN;
        S(S < -1000 | S < 0 | S > 45) = NaN;
        la(la < -1000) = NaN;
        lo(lo < -1000) = NaN;

        n = min([numel(tdt), numel(z), numel(T), numel(S), numel(la), numel(lo)]);
        Time  = [Time;  tdt(1:n)];
        Depth = [Depth; z(1:n)];
        Temp  = [Temp;  T(1:n)];
        Salt  = [Salt;  S(1:n)];
        Lat   = [Lat;   la(1:n)];
        Lon   = [Lon;   lo(1:n)];
        FileID = [FileID; repmat(i, n, 1)];
    end

    good = ~isnat(Time) & isfinite(Depth) & isfinite(Temp) & isfinite(Salt) & Depth >= 0;
    G = table(Time(good), Depth(good), Temp(good), Salt(good), Lat(good), Lon(good), FileID(good), ...
        'VariableNames', {'Time','Depth','Temperature','Salinity','Latitude','Longitude','FileID'});
end

function tdt = decodeGliderTime(t)
    t = double(t(:));
    medt = median(t, "omitnan");
    if medt > 1e8
        tdt = datetime(t, "ConvertFrom", "posixtime");
    elseif medt > 5e5
        tdt = datetime(t, "ConvertFrom", "datenum");
    else
        tdt = datetime(1970,1,1) + days(t);
    end
end

function P = makeSSP(G, period, zProfile)
    idx = periodMask(G, period);
    assert(any(idx), "No observations found for period %s.", period.name);
    P = buildSSPFromRows(G, find(idx), zProfile);
    P.name = period.name;
end

function idx = periodMask(G, period)
    if period.name == "Before"
        idx = G.Time < period.t1;
    elseif period.name == "After"
        idx = G.Time >= period.t0;
    else
        idx = G.Time >= period.t0 & G.Time <= period.t1;
    end
end

function cBoot = bootstrapSSP(G, period, zProfile, nBoot)
    idxPeriod = periodMask(G, period);
    ids = unique(G.FileID(idxPeriod));
    assert(~isempty(ids), "No bootstrap segments found for period %s.", period.name);

    rowsById = cell(numel(ids),1);
    for i = 1:numel(ids)
        rowsById{i} = find(idxPeriod & G.FileID == ids(i));
    end

    cBoot = nan(numel(zProfile), nBoot);
    for b = 1:nBoot
        pick = randi(numel(ids), numel(ids), 1);
        rows = vertcat(rowsById{pick});
        P = buildSSPFromRows(G, rows, zProfile);
        cBoot(:,b) = P.c;
    end
end

function P = buildSSPFromRows(G, rows, zProfile)
    z = G.Depth(rows);
    T = G.Temperature(rows);
    S = G.Salinity(rows);
    lat0 = median(G.Latitude(rows), "omitnan");
    edges = [zProfile - mean(diff(zProfile))/2; zProfile(end) + mean(diff(zProfile))/2];
    bin = discretize(z, edges);
    nb = numel(zProfile);

    Tbin = accumarray(bin(~isnan(bin)), T(~isnan(bin)), [nb 1], ...
        @(x) median(x, "omitnan"), NaN);
    Sbin = accumarray(bin(~isnan(bin)), S(~isnan(bin)), [nb 1], ...
        @(x) median(x, "omitnan"), NaN);
    Nbin = accumarray(bin(~isnan(bin)), 1, [nb 1], @sum, 0);

    Tbin(Nbin < 3) = NaN;
    Sbin(Nbin < 3) = NaN;
    Tbin = fillmissing(Tbin, "linear", "EndValues", "nearest");
    Sbin = fillmissing(Sbin, "linear", "EndValues", "nearest");

    c = mackenzieSoundSpeed(Tbin, Sbin, zProfile);

    P.name = "";
    P.z = zProfile;
    P.T = Tbin;
    P.S = Sbin;
    P.c = smoothdata(c, "movmean", 5, "omitnan");
    P.nObs = numel(rows);
    P.lat = lat0;
end

function c = mackenzieSoundSpeed(T, S, z)
    % Mackenzie (1981), T degC, S PSU, z m.
    c = 1448.96 ...
        + 4.591.*T ...
        - 5.304e-2.*T.^2 ...
        + 2.374e-4.*T.^3 ...
        + 1.340.*(S - 35) ...
        + 1.630e-2.*z ...
        + 1.675e-7.*z.^2 ...
        - 1.025e-2.*T.*(S - 35) ...
        - 7.139e-13.*T.*z.^3;
end

function rays = traceRaysSnell(zSSP, cSSP, zSrc, anglesDeg, maxRangeM, bottomM)
    dr = 80;
    r = 0:dr:maxRangeM;
    cSrc = interp1(zSSP, cSSP, zSrc, "linear", "extrap");
    rays = struct("angle", {}, "r", {}, "z", {});

    for ia = 1:numel(anglesDeg)
        theta0 = deg2rad(anglesDeg(ia));
        p = cos(theta0) / cSrc;
        sgn = sign(theta0);
        if sgn == 0
            sgn = 1;
        end

        zr = nan(size(r));
        zr(1) = zSrc;
        for ir = 2:numel(r)
            cc = interp1(zSSP, cSSP, zr(ir-1), "linear", "extrap");

            % Turning point when Snell's law would make cos(theta) > 1.
            if p * cc >= 0.99995
                sgn = -sgn;
                pc = 0.995; % move the ray away from the caustic instead of drawing a flat numerical stall
            else
                pc = max(min(p * cc, 0.99995), 0.0001);
            end

            slope = sgn * sqrt(max(1 - pc.^2, 0)) / pc;
            znew = zr(ir-1) + slope * dr;

            if znew < 0
                znew = -znew;
                sgn = abs(sgn);
            elseif znew > bottomM
                znew = 2*bottomM - znew;
                sgn = -abs(sgn);
            end
            zr(ir) = znew;
        end

        rays(ia).angle = anglesDeg(ia);
        rays(ia).r = r;
        rays(ia).z = zr;
    end
end

function [TL, z, r] = runMMPEstylePE(zSSP, cSSP, zSrc, freqHz, rangeM, dz, bottomM)
    z = (0:dz:bottomM).';
    r = rangeM(:).';
    c = interp1(zSSP, cSSP, z, "linear", "extrap");
    c = fillmissing(c, "nearest");

    c0 = interp1(z, c, zSrc, "linear", "extrap");
    k0 = 2*pi*freqHz/c0;
    n2 = (c0 ./ c).^2;

    interior = 2:numel(z)-1;
    zi = z(interior);
    ni2 = n2(interior);
    Ni = numel(interior);

    e = ones(Ni,1);
    D2 = spdiags([e -2*e e], -1:1, Ni, Ni) / dz^2;
    I = speye(Ni);
    A = I - 1i * mean(diff(r)) / (4*k0) * D2;
    B = I + 1i * mean(diff(r)) / (4*k0) * D2;

    lambda = c0 / freqHz;
    sigma = max(4*dz, 2.5*lambda);
    u = exp(-((zi - zSrc)/sigma).^2) - exp(-((zi + zSrc)/sigma).^2);
    u = u / max(abs(u));

    sponge = ones(Ni,1);
    spongeStart = bottomM - max(60, 0.18*bottomM);
    js = zi > spongeStart;
    sponge(js) = exp(-0.035*((zi(js) - spongeStart)/(bottomM - spongeStart)).^2);

    phase = exp(1i * k0 * mean(diff(r)) / 4 * (ni2 - 1));

    U = zeros(numel(z), numel(r), "single");
    U(interior,1) = single(u);

    for ir = 2:numel(r)
        u = phase .* u;
        u = A \ (B * u);
        u = phase .* u;
        u = u .* sponge;
        U(interior,ir) = single(u);
    end

    amp = double(abs(U));
    ref = max(amp(:,1), [], "omitnan");
    amp = max(amp/ref, 1e-10);

    spreading = 10*log10(max(r, 1)); % cylindrical spreading term, range in m
    TL = -20*log10(amp) + spreading;
    TL(:,1) = TL(:,2);
    TL = smoothdata(TL, 2, "movmean", 5);
end

function map = natureBlueRed(n)
    x = [0 0.25 0.5 0.75 1];
    c = [ ...
        0.075 0.224 0.416
        0.230 0.468 0.711
        0.925 0.930 0.900
        0.852 0.453 0.333
        0.600 0.110 0.100];
    map = interp1(x, c, linspace(0,1,n), "pchip");
    map = max(min(map,1),0);
end

function map = natureThermal(n)
    x = [0 0.18 0.38 0.62 0.82 1];
    c = [ ...
        0.070 0.110 0.260
        0.110 0.310 0.520
        0.220 0.560 0.610
        0.910 0.740 0.350
        0.850 0.360 0.200
        0.430 0.090 0.120];
    map = interp1(x, c, linspace(0,1,n), "pchip");
    map = max(min(map,1),0);
end

function panelLabel(ax, str)
    text(ax, 0.015, 0.965, str, "Units", "normalized", ...
        "HorizontalAlignment", "left", "VerticalAlignment", "top", ...
        "FontWeight", "bold", "FontSize", 13, "Color", "k", ...
        "BackgroundColor", "w", "Margin", 1.5);
end

function s = dateSpanText(t0, t1)
    s = sprintf("%s-%s", string(t0, "MM/dd HH:mm"), string(t1, "MM/dd HH:mm"));
end

function lim = paddedLimits(x, pad)
    x = x(isfinite(x));
    if isempty(x)
        lim = [0 1];
        return;
    end
    lo = floor((min(x)-pad)/5)*5;
    hi = ceil((max(x)+pad)/5)*5;
    if lo == hi
        hi = lo + 5;
    end
    lim = [lo hi];
end