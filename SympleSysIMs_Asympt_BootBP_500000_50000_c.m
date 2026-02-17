close all
clear
clc

% -----------------------------------------------------------
%  Read data from Excel
% ------------------------------------------------------------
filename1 = 'Process_EMRALD_DebugLog_500000_Emanuele.xlsx';
sheet     = 1;

% x labels
xlabels = readcell(filename1, 'FileType','spreadsheet', ...
                   'Sheet', sheet, 'Range','C7:F7');

% x data
xorig = readmatrix(filename1, 'FileType','spreadsheet', ...
                   'Sheet', sheet, 'Range','C8:F500007');
xorig(isnan(xorig)) = 0;

% y label
tmp_ylabel = readcell(filename1, 'FileType','spreadsheet', ...
                      'Sheet', sheet, 'Range','G7');
ylabel_str = tmp_ylabel{1};

% y data
yorig = readmatrix(filename1, 'FileType','spreadsheet', ...
                   'Sheet', sheet, 'Range','G8:G500007');
yorig(isnan(yorig)) = 0;

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
%  Settings
% ------------------------------------------------------------
% Large and small sample sizes
Nbig   = [500,5000,10000,15000,25000,50000,75000,150000,250000,500000];
Nsmall = [500,1000,2000,5000,10000,15000,20000,25000,30000,50000];

NBoot  = 100;

%% Bostrap Reliability Point Estimate


Unrel_big  = BootUnrMetric(xorig, yorig, Nbig, @precalcsb, NBoot)
Unrel_small= BootUnrMetric(xorig, yorig, Nsmall, @precalcsb, NBoot)


figure(1); set(gcf,'WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Big N
nexttile(1); hold on
for j = 1:1
    boxplot(squeeze(Unrel_big(:,:,j))');  % NBoot x length(Nbig)
end
xlabel('N'); title('Up to N=500,000');
xticks(1:numel(Nbig)); xticklabels(string(Nbig));
yline(0.0128,'--r','LineWidth',3);
ylim([0 0.05])
xtickangle(45); set(gca,'FontSize',24); box on; hold off

% Small N
nexttile(2); hold on
for j = 1
    boxplot(squeeze(Unrel_small(:,:,j))');
end
xlabel('N'); title('Up to N=50,000');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
yline(0.0128,'--r','LineWidth',3);
ylim([0 0.05])
xtickangle(45); set(gca,'FontSize',24); box on; hold off


filename = matlab.desktop.editor.getActiveFilename;
[~, filename, ~] = fileparts(filename)
saveas(gcf, ['../Figures/' filename '_Unrel.png']);
saveas(gcf, ['../Figures/' filename '_Unrel.epsc']);
saveas(gcf, ['../Figures/' filename '_Unrel.fig']);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importance Measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -----------------------------------------------------------
%  Bootstrap for all indices and both N-sets
% ------------------------------------------------------------
% Big N
FVb_big       = computeBootMetric(xorig, yorig, Nbig,   @FVf,       NBoot);
RAWb_big      = computeBootMetric(xorig, yorig, Nbig,   @RAWf,      NBoot);
RRWb_big      = computeBootMetric(xorig, yorig, Nbig,   @RRWf,      NBoot);
Birnbaumb_big = computeBootMetric(xorig, yorig, Nbig,   @Birnbaumf, NBoot);
BPb_big       = computeBootMetric(xorig, yorig, Nbig,   @BPf,       NBoot);

% Small N
FVb_small       = computeBootMetric(xorig, yorig, Nsmall, @FVf,       NBoot);
RAWb_small      = computeBootMetric(xorig, yorig, Nsmall, @RAWf,      NBoot);
RRWb_small      = computeBootMetric(xorig, yorig, Nsmall, @RRWf,      NBoot);
Birnbaumb_small = computeBootMetric(xorig, yorig, Nsmall, @Birnbaumf, NBoot);
BPb_small       = computeBootMetric(xorig, yorig, Nsmall, @BPf,       NBoot);

%% -----------------------------------------------------------
%  Means (last N for each set)
% ------------------------------------------------------------
mFVb_big       = squeeze(mean(FVb_big(end,:,:),       2));
mRAWb_big      = squeeze(mean(RAWb_big(end,:,:),      2));
mRRWb_big      = squeeze(mean(RRWb_big(end,:,:),      2));
mBirnbaumb_big = squeeze(mean(Birnbaumb_big(end,:,:), 2));
mBPb_big       = squeeze(mean(BPb_big(end,:,:),       2));

mFVb_small       = squeeze(mean(FVb_small(end,:,:),       2));
mRAWb_small      = squeeze(mean(RAWb_small(end,:,:),      2));
mRRWb_small      = squeeze(mean(RRWb_small(end,:,:),      2));
mBirnbaumb_small = squeeze(mean(Birnbaumb_small(end,:,:), 2));
mBPb_small       = squeeze(mean(BPb_small(end,:,:),       2));



% -----------------------------------------------------------
%  Boxplot Figures (big vs small N, 1x2 tiledlayout)
% ------------------------------------------------------------

% ---------- Fussell-Vesely ----------
figure(2); set(gcf,'WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Big N
nexttile(1); hold on
for j = 1:nComp
    boxplot(squeeze(FVb_big(:,:,j))');  % NBoot x length(Nbig)
end
xlabel('N'); title('Fussell-Vesely (big N)');
xticks(1:numel(Nbig)); xticklabels(string(Nbig));
yline(0.61,'--r','LineWidth',3);
yline(0.41,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

% Small N
nexttile(2); hold on
for j = 1:nComp
    boxplot(squeeze(FVb_small(:,:,j))');
end
xlabel('N'); title('Fussell-Vesely (small N)');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
yline(0.61,'--r','LineWidth',3);
yline(0.41,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

%% ---------- RAW ----------
figure(3); set(gcf,'WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile(1); hold on
for j = 1:nComp
    boxplot(squeeze(RAWb_big(:,:,j))');
end
xlabel('N'); title('RAW (big N)');
xticks(1:numel(Nbig)); xticklabels(string(Nbig));
yline(8.84,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

nexttile(2); hold on
for j = 1:nComp
    boxplot(squeeze(RAWb_small(:,:,j))');
end
xlabel('N'); title('RAW (small N)');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
yline(8.84,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

%% ---------- RRW ----------
figure(4); set(gcf,'WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile(1); hold on
for j = 1:nComp
    boxplot(squeeze(RRWb_big(:,:,j))');
end
xlabel('N'); title('RRW (big N)');
xticks(1:numel(Nbig)); xticklabels(string(Nbig));
xtickangle(45); set(gca,'FontSize',24); box on; hold off

nexttile(2); hold on
for j = 1:nComp
    boxplot(squeeze(RRWb_small(:,:,j))');
end
xlabel('N'); title('RRW (small N)');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
xtickangle(45); set(gca,'FontSize',24); box on; hold off

%% ---------- Birnbaum ----------
figure(5); set(gcf,'WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile(1); hold on
for j = 1:nComp
    boxplot(squeeze(Birnbaumb_big(:,:,j))');
end
xlabel('N'); title('Birnbaum (big N)');
ylim([0 0.4])
xticks(1:numel(Nbig)); xticklabels(string(Nbig));
yline(0.108,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

nexttile(2); hold on
for j = 1:nComp
    boxplot(squeeze(Birnbaumb_small(:,:,j))');
end
xlabel('N'); title('Birnbaum (small N)');
ylim([0 0.4])
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
yline(0.108,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

%% ---------- Barlow-Proschan ----------
figure(6); set(gcf,'WindowState','maximized');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile(1); hold on
for j = 1:nComp
    boxplot(squeeze(BPb_big(:,:,j))');
end
ylim([0 0.8]);
xlabel('N'); title('Barlow-Proschan (big N)');
xticks(1:numel(Nbig)); xticklabels(string(Nbig));
yline(0.3,'--r','LineWidth',3);
yline(0.2,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off

nexttile(2); hold on
for j = 1:nComp
    boxplot(squeeze(BPb_small(:,:,j))');
end
ylim([0 0.8]);
xlabel('N'); title('Barlow-Proschan (small N)');
xticks(1:numel(Nsmall)); xticklabels(string(Nsmall));
yline(0.3,'--r','LineWidth',3);
yline(0.2,'--r','LineWidth',3);
xtickangle(45); set(gca,'FontSize',24); box on; hold off



%% -----------------------------------------------------------
%  Overall horizontal-bar summary (small N set, as in your end)
% ------------------------------------------------------------
figure('Units','normalized','Position',[0.05 0.05 0.9 0.9]);

xlabels2  = arrayfun(@(i) sprintf('$C_{%d}$', i), 1:numel(xlabels), 'UniformOutput', false);
idxLabels = [1 2 3 4];   % choose which components to show as labels

metricsMeans = {mBirnbaumb_small, mRRWb_small, mRAWb_small, ...
                mFVb_small,     mBPb_small};
titlesCell   = {'Birnbaum','RRW','RAW','Fussell-Vesely','Barlow-Proschan'};

for i = 1:5
    subplot(2,3,i);
    barh(metricsMeans{i}(:), 'BarWidth', 0.85);
    if i == 5
        xlim([0 0.06]);
    end
    yticks(idxLabels);
    yticklabels(xlabels2(idxLabels));
    set(gca,'TickLabelInterpreter','latex','FontSize',22);
    title(titlesCell{i},'FontSize',24);
    box on;
end

%% ========================
%  Saving Values
%  =========================
save('SympleSysIMs_Asympt_BootBP_500000_50000_c.mat', 'Birnbaumb_big','Birnbaumb_small','BPb_big','BPb_small','FVb_big','FVb_small','RAWb_big','RAWb_small','RRWb_big','RRWb_small','xlabels','mBirnbaumb_big','mBirnbaumb_small','mBPb_big','mBPb_small','mFVb_big','mFVb_small','mRAWb_big','mRAWb_small','mRRWb_small','mRRWb_big');

%% ============================================================
%  Helper functions
% ============================================================


function UB = BootUnrMetric(xorig,yorig, Nvec, funHandle, NBoot)
    % Returns array: [numN x NBoot x nComp]
    nN    = numel(Nvec);
    nComp = 1;
    B     = zeros(nN, NBoot, nComp);

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

%% --------- Your original index functions --------------------


function [U]=precalcsb(x,y)
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
n=sum(ybin);
ni=sum(xbin);
U=n/size(y,1);
end


function [xbin,ybin,n,ni]=precalcs(x,y)
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
n=sum(ybin);
ni=sum(xbin);
end

function [FV]=FVf(x,y)
[xbin,ybin,n,ni]=precalcs(x,y);
xbinFV=zeros(size(x));
for j=1:size(x,2)
    indxjFV=find(and(y>0,x(:,j)>0));
    xbinFV(indxjFV,j)=1;
end
fi=sum(xbinFV); %Number of times the system fails given that component i has failed
FV=fi./n;
end

function [RAW]=RAWf(x,y)
[xbin,ybin,n,ni]=precalcs(x,y);
N=size(x,1);
for j=1:size(x,2)
    rawnum(j)=sum(ybin(find(xbin(:,j)==1)));
end
RAW=(rawnum./ni)/(n/N);
end

function [RRW]=RRWf(x,y)
[xbin,ybin,n,ni]=precalcs(x,y);
N=size(x,1);
for j=1:size(x,2)
    rrwnum(j)=sum(ybin(find(xbin(:,j)==0)));
end
RRW=(n/N)./(rrwnum./(N-ni));
end

function [Birnbaum]=Birnbaumf(x,y)
[xbin,ybin,n,ni]=precalcs(x,y);
N=size(x,1);
for j=1:size(x,2)
    rawnum(j)=sum(ybin(find(xbin(:,j)==1)));
    rrwnum(j)=sum(ybin(find(xbin(:,j)==0)));
end
Birnbaum=rawnum./ni-rrwnum./(N-ni);
end

function [BP]=BPf(x,y)
[xbin,ybin,n,ni]=precalcs(x,y);
xbinBP=zeros(size(xbin));
NBp=size(x,1);
for j=1:size(x,2)
    indxjBP=find(and(y>0,y==x(:,j)));
    xbinBP(indxjBP,j)=1;
    fbp(j)=sum(ybin(find(xbinBP(:,j)==1)));
end
BP=fbp/sum(ybin);
end
