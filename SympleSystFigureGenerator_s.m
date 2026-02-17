close all
clear
clc

%% -----------------------------------------------------------
%  Load previously saved workspace (all variables except x,y,xorig,yorig)
% ------------------------------------------------------------

dataFile = 'SympleSysIMs_Asympt_BootBP_500000_50000_DIM.mat';  
load(dataFile);
filename = matlab.desktop.editor.getActiveFilename;
[~, filename, ~] = fileparts(filename)
% Number of components (should be consistent with saved variables)
nComp = numel(xlabels);

% Script name for figure saving (this script's name)
nameforreproducibility1 = matlab.desktop.editor.getActiveFilename;
[~, nameforreproducibility, ~] = fileparts(nameforreproducibility1);

Nbig   = [500,5000,10000,15000,25000,50000,75000,150000,250000,500000];
Nsmall = [500,1000,2000,5000,10000,15000,20000,25000,30000,50000];

%% -----------------------------------------------------------
%  Boxplot Figures (big vs N up to 50000, 1x2 tiledlayout)
% ------------------------------------------------------------

% ---------- Birnbaum ----------
figure(1); set(gcf,'WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>

nexttile(1); hold on
for j = 1:nComp
    boxplot(squeeze(Birnbaumb_big(:,:,j))');
end
xlabel('N'); title('Birnbaum (N up to 500000)');
xticks(1:numel(Nbig)); xticklabels(string(Nbig));
ylim([0 0.5]);
yline(0.108,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

nexttile(2); hold on
for j = 1:nComp
    boxplot(squeeze(Birnbaumb_small(:,:,j))');
end
xlabel('N'); title('Birnbaum (N up to 50000)');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
yline(0.108,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off
saveas(gcf, ['../Figures/' filename '_Birnbaum.png']);
saveas(gcf, ['../Figures/' filename '_Birnbaum.epsc']);
saveas(gcf, ['../Figures/' filename '_Birnbaum.fig']);


% ---------- Barlow-Proschan ----------
figure(2); set(gcf,'WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>

nexttile(1); hold on
for j = 1:nComp
    boxplot(squeeze(BPb_big(:,:,j))');
end
ylim([0 0.8]);
xlabel('N'); title('Barlow-Proschan (N up to 500000)');
xticks(1:numel(Nbig)); xticklabels(string(Nbig));
yline(0.3,'--r','LineWidth',3);
yline(0.2,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

nexttile(2); hold on
for j = 1:nComp
    boxplot(squeeze(BPb_small(:,:,j))');
end
ylim([0 0.8]);
xlabel('N'); title('Barlow-Proschan (N up to 50000)');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
yline(0.3,'--r','LineWidth',3);
yline(0.2,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

saveas(gcf, ['../Figures/' filename '_BP.png']);
saveas(gcf, ['../Figures/' filename '_BP.epsc']);
saveas(gcf, ['../Figures/' filename '_BP.fig']);




% ---------- Fussell-Vesely ----------
figure(3); set(gcf,'WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>

% N up to 500000
nexttile(1); hold on
for j = 1:nComp
    boxplot(squeeze(FVb_big(:,:,j))');  % NBoot x length(Nbig)
end
xlabel('N'); title('Fussell-Vesely (N up to 500000)');
xticks(1:numel(Nbig)); xticklabels(string(Nbig));
yline(0.61,'--r','LineWidth',3);
yline(0.41,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

% N up to 50000
nexttile(2); hold on
for j = 1:nComp
    boxplot(squeeze(FVb_small(:,:,j))');
end
xlabel('N'); title('Fussell-Vesely (N up to 50000)');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
yline(0.61,'--r','LineWidth',3);
yline(0.41,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

saveas(gcf, ['../Figures/' filename '_FV.png']);
saveas(gcf, ['../Figures/' filename '_FV.epsc']);
saveas(gcf, ['../Figures/' filename '_FV.fig']);

% ---------- RAW ----------
figure(4); set(gcf,'WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>

nexttile(1); hold on
for j = 1:nComp
    boxplot(squeeze(RAWb_big(:,:,j))');
end
xlabel('N'); title('RAW (N up to 500000)');
xticks(1:numel(Nbig)); xticklabels(string(Nbig));
yline(8.84,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

nexttile(2); hold on
for j = 1:nComp
    boxplot(squeeze(RAWb_small(:,:,j))');
end
xlabel('N'); title('RAW (N up to 50000)');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
yline(8.84,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

saveas(gcf, ['../Figures/' filename '_RAW.png']);
saveas(gcf, ['../Figures/' filename '_RAW.epsc']);
saveas(gcf, ['../Figures/' filename '_RAW.fig']);


% ---------- RRW ----------
figure(5); set(gcf,'WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>

nexttile(1); hold on
for j = 1:nComp
    boxplot(squeeze(RRWb_big(:,:,j))');
end
yline(2.41,'--r','LineWidth',3);
yline(1.63,'--r','LineWidth',3);
xlabel('N'); title('RRW (N up to 500000)');
xticks(1:numel(Nbig)); xticklabels(string(Nbig));
xtickangle(45); set(gca,'FontSize',24); box on; hold off

nexttile(2); hold on
for j = 1:nComp
    boxplot(squeeze(RRWb_small(:,:,j))');
end
xlabel('N'); title('RRW (N up to 50000)');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
xtickangle(45); set(gca,'FontSize',24); box on; hold off
yline(2.41,'--r','LineWidth',3);
yline(1.63,'--r','LineWidth',3);
saveas(gcf, ['../Figures/' filename '_RRW.png']);
saveas(gcf, ['../Figures/' filename '_RRW.epsc']);
saveas(gcf, ['../Figures/' filename '_RRW.fig']);

%% ---------- DIM ----------
figure(6); set(gcf,'WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>

nexttile(1); hold on
for j = 1:nComp
    boxplot(squeeze(DIMb_big(:,:,j))');
end
yline(0.30,'--r','LineWidth',3);
yline(0.199,'--r','LineWidth',3);
xlabel('N'); title('DIM (N up to 500000)');
xticks(1:numel(Nbig)); xticklabels(string(Nbig));
xtickangle(45); set(gca,'FontSize',24); box on; hold off

nexttile(2); hold on
for j = 1:nComp
    boxplot(squeeze(DIMb_small(:,:,j))');
end
xlabel('N'); title('DIM (N up to 50000)');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
xtickangle(45); set(gca,'FontSize',24); box on; hold off
yline(0.30,'--r','LineWidth',3);
yline(0.199,'--r','LineWidth',3);
saveas(gcf, ['../Figures/' filename '_DIM.png']);
saveas(gcf, ['../Figures/' filename '_DIM.epsc']);
saveas(gcf, ['../Figures/' filename '_DIM.fig']);

%% -----------------------------------------------------------
%  Overall horizontal-bar summary (N up to 500000 set)
% ------------------------------------------------------------
figure('Units','normalized','Position',[0.05 0.05 0.9 0.9]);

xlabels2  = arrayfun(@(i) sprintf('$C_{%d}$', i), 1:numel(xlabels), 'UniformOutput', false);
idxLabels = [1 2 3 4];   % choose which components to show as labels

metricsMeans = {mBirnbaumb_big, mBPb_big, mFVb_big, mRAWb_big, mRRWb_big};
titlesCell   = {'Birnbaum','Barlow-Proschan','Fussell-Vesely','RAW','RRW'};

for i = 1:5
    subplot(2,3,i);
    barh(metricsMeans{i}(:), 'BarWidth', 0.85);
    %if i == 5
    %    xlim([0 0.06]);
    %end
    yticks(idxLabels);
    yticklabels(xlabels2(idxLabels));
    set(gca,'TickLabelInterpreter','latex','FontSize',22);
    title(titlesCell{i},'FontSize',24);
    box on;
end


set(gcf,'WindowState','maximized');
saveas(gcf, ['../Figures/' filename '_SmAll.png']);


%% -----------------------------------------------------------
%  Overall horizontal-bar summary (N up to 50000 set)
% ------------------------------------------------------------
figure('Units','normalized','Position',[0.05 0.05 0.9 0.9]);

xlabels2  = arrayfun(@(i) sprintf('$C_{%d}$', i), 1:numel(xlabels), 'UniformOutput', false);
idxLabels = [1 2 3 4];   % choose which components to show as labels

metricsMeans = {mBirnbaumb_small, mBPb_small, mFVb_small, mRAWb_small, mRRWb_small};
titlesCell   = {'Birnbaum','Barlow-Proschan','Fussell-Vesely','RAW','RRW'};

for i = 1:5
    subplot(2,3,i);
    barh(metricsMeans{i}(:), 'BarWidth', 0.85);
    %if i == 5
    %    xlim([0 0.06]);
    %end
    yticks(idxLabels);
    yticklabels(xlabels2(idxLabels));
    set(gca,'TickLabelInterpreter','latex','FontSize',22);
    title(titlesCell{i},'FontSize',24);
    box on;
end


set(gcf,'WindowState','maximized');
saveas(gcf, ['../Figures/' filename '_SmAll.png']);



%% -----------------------------------------------------------
%  Save all figures to ../Figures (same logic as original script)
% ------------------------------------------------------------
outputDir = fullfile('..', 'Figures');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

figList = findobj('Type', 'figure');
for i = 1:length(figList)
    fig = figList(i);
    figName = sprintf('%s_Figure_%d', ...
        nameforreproducibility, fig.Number);
    saveas(fig, fullfile(outputDir, [figName, '.png']));
end
