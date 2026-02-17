close all
clear
clc

% -----------------------------------------------------------
%  Read data from Excel
% ------------------------------------------------------------

filename1 = 'Process_EMRALD_Complex_Sys_50000_IMS.xlsx';

sheet = 1;  % change if your data is on a specific sheet name/number

% 1) Import the block C18:S50018 (we'll read header row + data row separately)
% 2) First row (C18:S18) -> xlabels
xlabels = readcell(filename1, 'FileType','spreadsheet', 'Sheet', sheet, 'Range','C18:S18');

% 3) Rest of the matrix (C19:S50018) -> x (numeric)
xorig = readmatrix(filename1, 'FileType','spreadsheet', 'Sheet', sheet, 'Range','C19:S50018');

% 4) Read T18:T50018
% 5) First cell (T18) -> ylabel
tmp_ylabel = readcell(filename1, 'FileType','spreadsheet', 'Sheet', sheet, 'Range','T18');
ylabel = tmp_ylabel{1};   % extract the scalar from the cell

% 6) Remainder (T19:T50018) -> y (numeric vector)
yorig = readmatrix(filename1, 'FileType','spreadsheet', 'Sheet', sheet, 'Range','T19:T50018');

% Sanity checks
assert(numel(xlabels) == size(xorig,2), ...
    'xlabels count (%d) must match x columns (%d).', ...
    numel(xlabels), size(xorig,2));
assert(numel(yorig) == size(xorig,1), ...
    'Length of y (%d) must match number of rows in x (%d).', ...
    numel(yorig), size(xorig,1));

% Script name for figure saving
nameforreproducibility1 = matlab.desktop.editor.getActiveFilename;
[~, nameforreproducibility, ~] = fileparts(nameforreproducibility1);

nComp = size(xorig,2);   % number of components

%% -----------------------------------------------------------
%  Settings (only small N available)
% ------------------------------------------------------------
Nsmall = [500,1000,2000,5000,10000,15000,20000,25000,30000,50000];
NBoot  = 100;

%% Bootstrap Reliability Point Estimate (only small N)
Unrel_small = BootUnrMetric(xorig, yorig, Nsmall, @precalcsb, NBoot);

% -----------------------------------------------------------
%  Reliability boxplot (only small N, single axes)
% -----------------------------------------------------------
%
figure(1); set(gcf,'WindowState','maximized');
hold on
% Only first (and only) component of Unrel_small
boxplot(squeeze(Unrel_small(:,:,1))');   % NBoot x length(Nsmall)
xlabel('N');
title('Unreliability, F');
xticks(1:numel(Nsmall));
xticklabels(string(Nsmall));
%yline(0.0128,'--r','LineWidth',3);
ylim([0.1 0.2])
xtickangle(45);
set(gca,'FontSize',24);
box on; hold off

% Save figure 1
filename = matlab.desktop.editor.getActiveFilename;
[~, filename, ~] = fileparts(filename);
saveas(gcf, ['../Figures/' filename '_Unrel.png']);
saveas(gcf, ['../Figures/' filename '_Unrel.epsc']);
saveas(gcf, ['../Figures/' filename '_Unrel.fig']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importance Measures (only small N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -----------------------------------------------------------
%  Bootstrap for all indices (small N only)
% ------------------------------------------------------------
FVb_small       = computeBootMetric(xorig, yorig, Nsmall, @FVf,       NBoot);
RAWb_small      = computeBootMetric(xorig, yorig, Nsmall, @RAWf,      NBoot);
RRWb_small      = computeBootMetric(xorig, yorig, Nsmall, @RRWf,      NBoot);
Birnbaumb_small = computeBootMetric(xorig, yorig, Nsmall, @Birnbaumf, NBoot);
BPb_small       = computeBootMetric(xorig, yorig, Nsmall, @BPf,       NBoot);

%% -----------------------------------------------------------
%  Means (last N of small set)
% ------------------------------------------------------------
mFVb_small       = squeeze(mean(FVb_small(end,:,:),       2));
mRAWb_small      = squeeze(mean(RAWb_small(end,:,:),      2));
mRRWb_small      = squeeze(mean(RRWb_small(end,:,:),      2));
mBirnbaumb_small = squeeze(mean(Birnbaumb_small(end,:,:), 2));
mBPb_small       = squeeze(mean(BPb_small(end,:,:),       2));

% -----------------------------------------------------------
%  Boxplot Figures (small N only, one panel per figure)
% ------------------------------------------------------------

% ---------- Fussell-Vesely ----------
figure(2); set(gcf,'WindowState','maximized');
hold on
for j = 1:nComp
    boxplot(squeeze(FVb_small(:,:,j))');  % NBoot x length(Nsmall)
end
xlabel('N');
title('Fussell-Vesely');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
ylim([0 0.6])
%yline(0.61,'--r','LineWidth',3);
%yline(0.41,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

%% ---------- RAW ----------
figure(3); set(gcf,'WindowState','maximized');
hold on
for j = 1:nComp
    boxplot(squeeze(RAWb_small(:,:,j))');
end
xlabel('N'); title('RAW');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
%yline(8.84,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

%% ---------- RRW ----------
figure(4); set(gcf,'WindowState','maximized');
hold on
for j = 1:nComp
    boxplot(squeeze(RRWb_small(:,:,j))');
end
xlabel('N'); title('RRW');
ylim([0.7 2])
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
xtickangle(45); set(gca,'FontSize',24); box on; hold off

%% ---------- Birnbaum ----------
figure(5); set(gcf,'WindowState','maximized');
hold on
for j = 1:nComp
    boxplot(squeeze(Birnbaumb_small(:,:,j))');
end
xlabel('N'); title('Birnbaum');
ylim([0 0.8])
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
%yline(0.108,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

%% ---------- Barlow-Proschan ----------
figure(6); set(gcf,'WindowState','maximized');
hold on
for j = 1:nComp
    boxplot(squeeze(BPb_small(:,:,j))');
end
ylim([0 0.6]);
xlabel('N'); title('Barlow-Proschan');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
%yline(0.3,'--r','LineWidth',3);
%yline(0.2,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

%% -----------------------------------------------------------
%  Overall horizontal-bar summary (small N set)
% ------------------------------------------------------------
figure('Units','normalized','Position',[0.05 0.05 0.9 0.9]);

xlabels2  = arrayfun(@(i) sprintf('$C_{%d}$', i), 1:numel(xlabels), 'UniformOutput', false);
%idxLabels = [1 2 3 4];   % choose which components to show as labels

metricsMeans = {mBirnbaumb_small, mRRWb_small, mRAWb_small, ...
                mFVb_small,     mBPb_small};
titlesCell   = {'Birnbaum','RRW','RAW','Fussell-Vesely','Barlow-Proschan'};

for i = 1:5
    subplot(2,3,i);
    barh(metricsMeans{i}(:), 'BarWidth', 0.85);
    yticks([1:17])
    yticklabels(xlabels);
    title(titlesCell{i},'FontSize',24);
    box on;
end

%% ========================
%  Saving Values (small N only)
%  =========================
save(filename, ...
    'Birnbaumb_small','BPb_small','FVb_small','RAWb_small','RRWb_small', ...
    'xlabels', ...
    'mBirnbaumb_small','mBPb_small','mFVb_small','mRAWb_small','mRRWb_small');

%% ============================================================
%  Helper functions
% ============================================================

function UB = BootUnrMetric(xorig,yorig, Nvec, funHandle, NBoot)
    % Returns array: [numN x NBoot x nComp]
    nN    = numel(Nvec);
    nComp = 1;  % precalcsb returns a scalar
    UB    = zeros(nN, NBoot, nComp);

    for u = 1:nN
        x = xorig(1:Nvec(u),:);
        y = yorig(1:Nvec(u));
        bootfun = @(x,y) funHandle(x,y);
        tmp = bootstrp(NBoot, bootfun, x, y);  % NBoot x nComp
        UB(u,:,:) = tmp;
    end
end

function B = computeBootMetric(xorig, yorig, Nvec, funHandle, NBoot)
    % Returns array: [numN x NBoot x nComp]
    nN    = numel(Nvec);
    nComp = size(xorig,2);
    B     = zeros(nN, NBoot, nComp);

    for u = 1:nN
        x = xorig(1:Nvec(u),:);
        y = yorig(1:Nvec(u));
        bootfun = @(x,y) funHandle(x,y);
        tmp = bootstrp(NBoot, bootfun, x, y);  % NBoot x nComp
        B(u,:,:) = tmp;
    end
end

%% --------- Index functions --------------------

function [U] = precalcsb(x,y)
    xbin=zeros(size(x));
    for j=1:size(x,2)
        indxj=find(x(:,j)>0);
        xbin(indxj,j)=1;
    end
    ybin=zeros(size(y));
    for j=1:size(y,2)
        indyj=find(y(:,j)>0);
        ybin(indyj,j)=1;
    end
    n = sum(ybin);
    ni = sum(xbin); %#ok<NASGU>
    U = n/size(y,1);
end

function [xbin,ybin,n,ni] = precalcs(x,y)
    xbin=zeros(size(x));
    for j=1:size(x,2)
        indxj=find(x(:,j)>0);
        xbin(indxj,j)=1;
    end
    ybin=zeros(size(y));
    for j=1:size(y,2)
        indyj=find(y(:,j)>0);
        ybin(indyj,j)=1;
    end
    n  = sum(ybin);
    ni = sum(xbin);
end

function [FV] = FVf(x,y)
    [xbin,ybin,n,ni] = precalcs(x,y);
    xbinFV=zeros(size(x));
    for j=1:size(x,2)
        indxjFV=find(and(y>0,x(:,j)>0));
        xbinFV(indxjFV,j)=1;
    end
    fi = sum(xbinFV); %Number of times the system fails given that component i has failed
    FV = fi./n;
end

function [RAW] = RAWf(x,y)
    [xbin,ybin,n,ni] = precalcs(x,y);
    N = size(x,1);
    rawnum = zeros(1,size(x,2));
    for j=1:size(x,2)
        rawnum(j) = sum(ybin(find(xbin(:,j)==1)));
    end
    RAW = (rawnum./ni)/(n/N);
end

function [RRW] = RRWf(x,y)
    [xbin,ybin,n,ni] = precalcs(x,y);
    N = size(x,1);
    rrwnum = zeros(1,size(x,2));
    for j=1:size(x,2)
        rrwnum(j) = sum(ybin(find(xbin(:,j)==0)));
    end
    RRW = (n/N)./(rrwnum./(N-ni));
end

function [Birnbaum] = Birnbaumf(x,y)
    [xbin,ybin,n,ni] = precalcs(x,y);
    N = size(x,1);
    rawnum = zeros(1,size(x,2));
    rrwnum = zeros(1,size(x,2));
    for j=1:size(x,2)
        rawnum(j) = sum(ybin(find(xbin(:,j)==1)));
        rrwnum(j) = sum(ybin(find(xbin(:,j)==0)));
    end
    Birnbaum = rawnum./ni - rrwnum./(N-ni);
end

function [BP] = BPf(x,y)
    [xbin,ybin,n,ni] = precalcs(x,y); %#ok<ASGLU>
    xbinBP=zeros(size(xbin));
    NBp = size(x,1); %#ok<NASGU>
    fbp = zeros(1,size(x,2));
    for j=1:size(x,2)
        indxjBP=find(and(y>0,y==x(:,j)));
        xbinBP(indxjBP,j)=1;
        fbp(j)=sum(ybin(find(xbinBP(:,j)==1)));
    end
    BP = fbp/sum(ybin);
end
