close all
clear
clc

%% -----------------------------------------------------------
%  Load previously saved workspace (all variables except x,y,xorig,yorig)
% ------------------------------------------------------------
dataFile = 'ComplexSysIMs_Asympt_Boot_50000_only.mat';  
load(dataFile);
filename = matlab.desktop.editor.getActiveFilename;
[~, filename, ~] = fileparts(filename)
% Number of components (should be consistent with saved variables)
nComp = numel(xlabels);

% Script name for figure saving (this script's name)
nameforreproducibility1 = matlab.desktop.editor.getActiveFilename;
[~, nameforreproducibility, ~] = fileparts(nameforreproducibility1);

Nsmall = [500,1000,2000,5000,10000,15000,20000,25000,30000,50000];
idxLabels = [1:17];
metricsMeans = {mBirnbaumb_small, mDIM_small, mBPb_small, mFVb_small, mRAWb_small, mRRWb_small};
titlesCell   = {'Birnbaum','Barlow-Proschan','Fussell-Vesely','RAW','RRW'};
Ntarget = 50000;

%% -----------------------------------------------------------
%  Boxplot Figures (only N up to 50000, NO big sample)
% ------------------------------------------------------------

% ---------- Birnbaum ----------
figure(1); 
subplot(1,2,1)
hold on
for j = 1:nComp
    boxplot(squeeze(Birnbaumb_small(:,:,j))');
end
xlabel('N'); title('$B_i$, Bootstrap Intervals','Interpreter','latex');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
xtickangle(45); set(gca,'FontSize',24); 
subplot(1,2,2)
    barh(metricsMeans{1}(:), 'BarWidth', 0.85);
    yticks(idxLabels);
    yticklabels(xlabels);
    set(gca,'TickLabelInterpreter','none');
    set(gca,'FontSize',18);
    box on;
title('$B_i$, Bootstrap Mean','FontSize',24,'Interpreter','latex');
set(gcf,'WindowState','maximized');

saveas(gcf,['../Figures/' filename '_Birnbaum.png']);

%saveas(gcf, ['../Figures/' filename '_Birnbaum.epsc']);
%saveas(gcf, ['../Figures/' filename '_Birnbaum.fig']);


%% ---------- Barlow-Proschan ----------
figure(2); set(gcf,'WindowState','maximized');

subplot(1,2,1)
hold on
for j = 1:nComp
    boxplot(squeeze(BPb_small(:,:,j))');
end
xlabel('N'); title('$BP_i$, Bootstrap Intervals','Interpreter','latex');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
ylim([0 0.6])
xtickangle(45); set(gca,'FontSize',24); 
subplot(1,2,2)
    barh(metricsMeans{2}(:), 'BarWidth', 0.85);
    yticks(idxLabels);
    yticklabels(xlabels);
    set(gca,'TickLabelInterpreter','none');
    set(gca,'FontSize',18);
    box on;
title('$BP_i$, Bootstrap Mean','FontSize',24,'Interpreter','latex');
set(gcf,'WindowState','maximized');

saveas(gcf,['../Figures/' filename '_BP.png']);

%saveas(gcf, ['../Figures/' filename '_BP.epsc']);
%saveas(gcf, ['../Figures/' filename '_BP.fig']);


%% ---------- Fussell-Vesely ----------
figure(3); set(gcf,'WindowState','maximized');

subplot(1,2,1)
hold on
for j = 1:nComp
    boxplot(squeeze(FVb_small(:,:,j))');
end
xlabel('N'); title('$FV_i$, Bootstrap Intervals','Interpreter','latex');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
ylim([0 0.5])
yline([0.005],'--r','LineWidth',3)
xtickangle(45); set(gca,'FontSize',24); 
subplot(1,2,2)
    barh(metricsMeans{3}(:), 'BarWidth', 0.85);
    yticks(idxLabels);
    yticklabels(xlabels);
    set(gca,'TickLabelInterpreter','none');
    set(gca,'FontSize',18);
    box on;
title('$FV_i$, Bootstrap Mean','FontSize',24,'Interpreter','latex');
set(gcf,'WindowState','maximized');

saveas(gcf,['../Figures/' filename '_FV.png']);

%saveas(gcf, ['../Figures/' filename '_BP.epsc']);
%saveas(gcf, ['../Figures/' filename '_BP.fig']);



%% ---------- RAW ----------
figure(4); set(gcf,'WindowState','maximized');

subplot(1,2,1)
hold on
for j = 1:nComp
    boxplot(squeeze(RAWb_small(:,:,j))');
end
xlabel('N'); title('$RAW_i$, Bootstrap Intervals','Interpreter','latex');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
ylim([0 max(max(max(RAWb_small)))])
yline(2,'--r','LineWidth',3)
xtickangle(45); set(gca,'FontSize',24); 
subplot(1,2,2)
    barh(metricsMeans{4}(:), 'BarWidth', 0.85);
    yticks(idxLabels);
    yticklabels(xlabels);
    set(gca,'TickLabelInterpreter','none');
    set(gca,'FontSize',18);
    box on;
title('$RAW_i$, Bootstrap Mean','FontSize',24,'Interpreter','latex');
set(gcf,'WindowState','maximized');

saveas(gcf,['../Figures/' filename '_RAW.png']);

%saveas(gcf, ['../Figures/' filename '_BP.epsc']);
%saveas(gcf, ['../Figures/' filename '_BP.fig']);


%% ---------- RRW ----------
figure(5); set(gcf,'WindowState','maximized');

subplot(1,2,1)
hold on
for j = 1:nComp
    boxplot(squeeze(RRWb_small(:,:,j))');
end
xlabel('N'); title('$RRW_i$, Bootstrap Intervals','Interpreter','latex');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
ylim([.85 max(max(max(RRWb_small)))])
xtickangle(45); set(gca,'FontSize',24); 
subplot(1,2,2)
    barh(metricsMeans{5}(:), 'BarWidth', 0.85);
    yticks(idxLabels);
    yticklabels(xlabels);
    set(gca,'TickLabelInterpreter','none');
    set(gca,'FontSize',18);
    box on;
title('$RRW_i$, Bootstrap Mean','FontSize',24,'Interpreter','latex');
set(gcf,'WindowState','maximized');

saveas(gcf,['../Figures/' filename '_RRW.png']);

%saveas(gcf, ['../Figures/' filename '_BP.epsc']);
%saveas(gcf, ['../Figures/' filename '_BP.fig']);

% ---------- Birnbaum ----------
figure(6); 
subplot(1,2,1)
hold on
for j = 1:nComp
    boxplot(squeeze(DIMb_small(:,:,j))');
end
xlabel('N'); title('$DIM_i$, Bootstrap Intervals','Interpreter','latex');
ylim([0 0.5])
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
xtickangle(45); set(gca,'FontSize',24); 
subplot(1,2,2)
    barh(metricsMeans{1}(:), 'BarWidth', 0.85);
    yticks(idxLabels);
    yticklabels(xlabels);
    set(gca,'TickLabelInterpreter','none');
    set(gca,'FontSize',18);
    box on;
title('$DIM_i$, Bootstrap Mean','FontSize',24,'Interpreter','latex');
set(gcf,'WindowState','maximized');

saveas(gcf,['../Figures/' filename '_DIM.png']);

%% -----------------------------------------------------------
%  Overall horizontal-bar summary (ONLY small sample)
% ------------------------------------------------------------
figure('Units','normalized','Position',[0.05 0.05 0.9 0.9]);

xlabels2  = arrayfun(@(i) sprintf('$C_{%d}$', i), 1:numel(xlabels), 'UniformOutput', false);
idxLabels = [1:17];

metricsMeans = {mBirnbaumb_small, mBPb_small, mDIM_small, mFVb_small, mRAWb_small, mRRWb_small};
titlesCell   = {'Birnbaum','Barlow-Proschan','DIM','Fussell-Vesely','RAW','RRW'};

for i = 1:6
    subplot(2,3,i);
    barh(metricsMeans{i}(:), 'BarWidth', 0.85);
    if or(i==1,i==4)
    yticks(idxLabels);
    yticklabels(xlabels);
    set(gca,'TickLabelInterpreter','none');
    end
    title(titlesCell{i},'FontSize',24);
    set(gca,'FontSize',16);
    box on;
end

set(gcf,'WindowState','maximized');
saveas(gcf, ['../Figures/' filename '_SmAll.png']);


%% Horizontal Pictures 



% Birnbaum
plot_bar_mean_ci_atN( ...
    Birnbaumb_small, Nsmall, xlabels, Ntarget, ...
    'Birnbaum', 'B_i', ...
    ['../Figures/' filename '_Birnbaum_barCI.png'] );

% Barlow Proschan
plot_bar_mean_ci_atN( ...
    BPb_small, Nsmall, xlabels, Ntarget, ...
    'Barlow-Proschan', 'BPi', ...
    ['../Figures/' filename '_BP_barCI.png'] );


% DIM
plot_bar_mean_ci_atN( ...
    DIMb_small, Nsmall, xlabels, Ntarget, ...
    'Differential Importance Measure', 'DIM_i', ...
    ['../Figures/' filename '_DIM_barCI.png'] );

% Fussell-Vesely
plot_bar_mean_ci_atN( ...
    FVb_small, Nsmall, xlabels, Ntarget, ...
    'Fussell-Vesely', 'FV_i', ...
    ['../Figures/' filename '_FV_barCI.png'] );

%RAWFussell-Vesely
plot_bar_mean_ci_atN( ...
    RAWb_small, Nsmall, xlabels, Ntarget, ...
    'Risk Achievement Worth', 'RAW_i', ...
    ['../Figures/' filename '_RAW_barCI.png'] );
% 
%RRW
plot_bar_mean_ci_atN( ...
    RRWb_small, Nsmall, xlabels, Ntarget, ...
    'Risk Reduction Worth', 'RRW_i', ...
    ['../Figures/' filename '_RRW_barCI.png'] );
% 

function plot_bar_mean_ci_atN(metric_small, Nsmall, xlabels, Ntarget, plotTitle, yLabel, outPng)
    % metric_small: [nN x nBoot x nComp]
    % xlabels: cell array of component names
    % Ntarget: e.g. 50000
    % outPng: e.g. '../Figures/myfile_Birnbaum.png'

    iN = find(Nsmall == Ntarget, 1, 'first');
    if isempty(iN)
        error('Ntarget=%g not found in Nsmall.', Ntarget);
    end

    % Extract bootstrap samples at Ntarget: [nBoot x nComp]
    S = squeeze(metric_small(iN,:,:));   % -> [nBoot x nComp]
    if isvector(S)
        S = S(:).'; % edge case: one component
    end

    m  = mean(S, 1, 'omitnan');
    ci = prctile(S, [2.5 97.5], 1);      % 2 x nComp
    lo = ci(1,:);
    hi = ci(2,:);

    nComp = numel(m);
    x = 1:nComp;

    figure('Color','w'); set(gcf,'WindowState','maximized');
    b = bar(x, m, 'BarWidth', 0.75);
    hold on;

    % Asymmetric error bars
    errLow  = m - lo;
    errHigh = hi - m;
    errorbar(x, m, errLow, errHigh, 'k', 'LineStyle','none', 'LineWidth', 1.8);

    grid on; box on;
    title(plotTitle, 'Interpreter','none'); % or 'latex' if you prefer
    if ~isempty(yLabel), ylabel(yLabel); end

    xticks(x);
    xticklabels(xlabels);
    set(gca,'TickLabelInterpreter','none');
    xtickangle(45);
    set(gca,'FontSize',22);

    % Ensure output folder exists
    outDir = fileparts(outPng);
    if ~isempty(outDir) && ~exist(outDir,'dir')
        mkdir(outDir);
    end

    exportgraphics(gcf, outPng, 'Resolution', 300);
end

%% Two metric plots

% B and BP
plot_grouped2_bar_mean_ci_atN( Birnbaumb_small, BPb_small, Nsmall, xlabels, Ntarget, ...
    'Birnbaum and Barlow Proschan', '$B_i$', '$BP_i$',...
    ['../Figures/' filename '_B_BP.png'] );

%% B and BP
plot_grouped2_bar_mean_ci_atN( DIMb_small, FVb_small, Nsmall, xlabels, Ntarget, ...
    'Dfferential Importance and Fussel-Vesely', '$DIM_i$', '$FV_i$',...
    ['../Figures/' filename '_DIM_FV.png'] );
%%
% RAW and RRW
plot_grouped2_bar_mean_ci_atN(RAWb_small, RRWb_small, Nsmall, xlabels, Ntarget, ...
    'RAW and RRW', '$RAW_i$', '$RRW_i$',...
    ['../Figures/' filename '_RAW_RRW.png'] );

function plot_grouped2_bar_mean_ci_atN(metricA_small, metricB_small, Nsmall, xlabels, Ntarget, ...
                                       titleStr, legendA, legendB, outPng)

    iN = find(Nsmall == Ntarget, 1, 'first');
    if isempty(iN)
        error('Ntarget=%g not found in Nsmall.', Ntarget);
    end

    SA = squeeze(metricA_small(iN,:,:)); % [nBoot x nComp]
    SB = squeeze(metricB_small(iN,:,:)); % [nBoot x nComp]

    mA  = mean(SA,1,'omitnan');  ciA = prctile(SA,[2.5 97.5],1);
    mB  = mean(SB,1,'omitnan');  ciB = prctile(SB,[2.5 97.5],1);

    loA = ciA(1,:); hiA = ciA(2,:);
    loB = ciB(1,:); hiB = ciB(2,:);

    nComp = numel(mA);
    x = 1:nComp;

    figure('Color','w'); set(gcf,'WindowState','maximized');

    Y = [mA(:) mB(:)];                % nComp x 2
    bh = bar(x, Y, 'grouped');        % two bars per component
    bh(1).FaceColor = [0.00 0.00 1.00];
    bh(2).FaceColor = [0.83 1.00 1.00];
    hold on;

    % Compute bar centers for errorbars
    xA = bh(1).XEndPoints;
    xB = bh(2).XEndPoints;

    errorbar(xA, mA, mA-loA, hiA-mA, 'k', 'LineStyle','none', 'LineWidth', 1.6);
    errorbar(xB, mB, mB-loB, hiB-mB, 'k', 'LineStyle','none', 'LineWidth', 1.6);

    grid on; box on;
    title(titleStr, 'Interpreter','none');
    ylabel({legendA, legendB},'Interpreter','latex')
    legend({legendA, legendB},'Interpreter','latex', ...
       'Location','northwest', ...
       'Box','on') %'Location','southoutside', 'Orientation','horizontal');

    xticks(x);
    xticklabels(xlabels);
    set(gca,'TickLabelInterpreter','none');
    xtickangle(45);
    set(gca,'FontSize',20);

    outDir = fileparts(outPng);
    if ~isempty(outDir) && ~exist(outDir,'dir')
        mkdir(outDir);
    end
    exportgraphics(gcf, outPng, 'Resolution', 300);
    ComplexSystFigureGenerator_B_BP
end

%% Three at a time
plot_grouped3_bar_mean_ci_atN( Birnbaumb_small, BPb_small, DIMb_small, Nsmall, xlabels, Ntarget, ...
    'Birnbaum, Barlow Proschan and DIM', 'B_i', 'BP_i','DIM_i',['../Figures/' filename '_B_BP_DIM.png'] );


function plot_grouped3_bar_mean_ci_atN(metricA_small, metricB_small, metricC_small, ...
                                       Nsmall, xlabels, Ntarget, ...
                                       titleStr, legendA, legendB, legendC, outPng)

    iN = find(Nsmall == Ntarget, 1, 'first');
    if isempty(iN)
        error('Ntarget=%g not found in Nsmall.', Ntarget);
    end

    % Extract bootstrap samples at Ntarget: [nBoot x nComp]
    SA = squeeze(metricA_small(iN,:,:));
    SB = squeeze(metricB_small(iN,:,:));
    SC = squeeze(metricC_small(iN,:,:));

    % Means and CIs
    mA  = mean(SA,1,'omitnan');  ciA = prctile(SA,[2.5 97.5],1);  loA = ciA(1,:); hiA = ciA(2,:);
    mB  = mean(SB,1,'omitnan');  ciB = prctile(SB,[2.5 97.5],1);  loB = ciB(1,:); hiB = ciB(2,:);
    mC  = mean(SC,1,'omitnan');  ciC = prctile(SC,[2.5 97.5],1);  loC = ciC(1,:); hiC = ciC(2,:);

    nComp = numel(mA);
    x = 1:nComp;

    figure('Color','w');
    set(gcf,'WindowState','maximized');

    % Grouped bars
    Y  = [mA(:) mB(:) mC(:)];     % nComp x 3
    bh = bar(x, Y, 'grouped');
    hold on;

    % Colors (match your style; adjust if you want)
    bh(1).FaceColor = [0 0 1];    % blue
    bh(2).FaceColor = [1 0 0];    % red
    bh(3).FaceColor = [0 0.6 0];  % green

    % Bar centers for errorbars
    xA = bh(1).XEndPoints;
    xB = bh(2).XEndPoints;
    xC = bh(3).XEndPoints;

    % Error bars (asymmetric)
    errorbar(xA, mA, mA-loA, hiA-mA, 'k', 'LineStyle','none', 'LineWidth', 1.6);
    errorbar(xB, mB, mB-loB, hiB-mB, 'k', 'LineStyle','none', 'LineWidth', 1.6);
    errorbar(xC, mC, mC-loC, hiC-mC, 'k', 'LineStyle','none', 'LineWidth', 1.6);

    grid on; box on;
    title(titleStr, 'Interpreter','none');

    legend({legendA, legendB, legendC}, ...
        'Interpreter','latex', ...
        'Location','northwest', ...
        'Box','on');

    xticks(x);
    xticklabels(xlabels);
    set(gca,'TickLabelInterpreter','none');
    xtickangle(45);
    set(gca,'FontSize',20);

    % Optional: reduce bottom whitespace (works well with rotated labels)
    ax = gca;
    ax.Position(2) = 0.18;
    ax.Position(4) = 0.72;
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))

    % Save
    outDir = fileparts(outPng);
    if ~isempty(outDir) && ~exist(outDir,'dir')
        mkdir(outDir);
    end
    exportgraphics(gcf, outPng, 'Resolution', 300);
    outPng
end
%% Safety Plane
%% ============================================================
%  Risk-informed plane (FV vs RAW) using MEAN bootstrap values
%  - Points = components
%  - Dashed red thresholds split into 4 RISC regions
%  - Labels + optional component names
% ============================================================



% ---- Mean bootstrap values (you already have these as *_small in your file) ----
FV  = mFVb_small(:);
RAW = mRAWb_small(:);

% ---- Component labels ----
if exist('xlabels','var') && numel(xlabels)==numel(FV)
    compNames = cellstr(string(xlabels));
else
    compNames = arrayfun(@(k) sprintf('Comp_%d',k), 1:numel(FV), 'UniformOutput', false);
end

nComp = numel(FV);

% ---- Thresholds (match your figure; change if you want different cutoffs) ----
xThr = 0.005;   % FV vertical threshold
yThr = 2.0;    % RAW horizontal threshold

% ---- Figure ----
figure('Color','w');
hold on; grid on; box on;

% Scatter points (diamond markers like your screenshot)
plot(FV, RAW, 'd', 'MarkerSize', 7,'LineWidth',5,'Color','r');

% Threshold lines (dashed red)
xline(xThr, '--r', 'LineWidth', 2);
yline(yThr, '--r', 'LineWidth', 2);

% Axis labels
xlabel('FV');
ylabel('RAW');

% Limits (auto but with some padding)
xmin = min(FV); xmax = max(FV);
ymin = min(RAW); ymax = max(RAW);
xpad = 0.05*(xmax - xmin + eps);
ypad = 0.05*(ymax - ymin + eps);
xlim([max(0,xmin-xpad), xmax+xpad]);
ylim([max(0,ymin-ypad), ymax+ypad]);

% ---- Region labels (placed relative to axes limits) ----
xl = xlim; yl = ylim;

% Quadrants based on thresholds:
% RISC 1: FV < xThr, RAW > yThr
% RISC 2: FV > xThr, RAW > yThr
% RISC 3: FV < xThr, RAW < yThr
% RISC 4: FV > xThr, RAW < yThr
text(xl(1)+0.02*(xl(2)-xl(1)), yl(1)+0.80*(yl(2)-yl(1)), 'RISC 1', ...
    'Color','b','FontWeight','bold','FontSize',24);
text(xl(1)+0.55*(xl(2)-xl(1)), yl(1)+0.78*(yl(2)-yl(1)), 'RISC 2', ...
    'Color','b','FontWeight','bold','FontSize',24);
text(xl(1)+0.02*(xl(2)-xl(1)), yl(1)+0.08*(yl(2)-yl(1)), 'RISC 3', ...
    'Color','b','FontWeight','bold','FontSize',24);
text(xl(1)+0.55*(xl(2)-xl(1)), yl(1)+0.08*(yl(2)-yl(1)), 'RISC 4', ...
    'Color','b','FontWeight','bold','FontSize',24);

title('Risk-informed plane using mean bootstrap values');

% ---- Optional: label each point with its component name (toggle on/off) ----
doLabelPoints = true;
if doLabelPoints
    for i = 1:nComp
        % small offset so text doesn't sit on marker
        text(FV(i) + 0.004*(xl(2)-xl(1)), RAW(i) + 0.01*(yl(2)-yl(1)), compNames{i}, ...
            'FontSize',20, 'Interpreter','none');
    end
end

% ---- Optional: list components by region in command window ----
idxR1 = (FV <= xThr) & (RAW >  yThr);
idxR2 = (FV >  xThr) & (RAW >  yThr);
idxR3 = (FV <= xThr) & (RAW <= yThr);
idxR4 = (FV >  xThr) & (RAW <= yThr);

fprintf('\n--- Components by region (FV threshold = %.3g, RAW threshold = %.3g) ---\n', xThr, yThr);
fprintf('RISC 1 (low FV, high RAW):\n'); disp(compNames(idxR1));
fprintf('RISC 2 (high FV, high RAW):\n'); disp(compNames(idxR2));
fprintf('RISC 3 (low FV, low RAW):\n'); disp(compNames(idxR3));
fprintf('RISC 4 (high FV, low RAW):\n'); disp(compNames(idxR4));

hold off;
set(gca,'FontSize',24)
set(gcf,'WindowState','maximized')
exportgraphics(gcf, ['../Figures/' filename 'RISCPlane.png'], 'Resolution', 300);
%['../Figures/' filename '_B_BP_DIM.png']