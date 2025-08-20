%% Plot sstructural stats

start_trees;
addpath('/cubric/data/c1441567/LargeScaleAnalysis2/complete_summary_stats')
load('/cubric/data/c1441567/LargeScaleAnalysis2/complete_summary_stats/total_stats.mat')

celltypes = {'microglia' 'astrocyte' 'oligodendrocyte' 'pyramidal' 'granule' 'purkinje' 'basket' 'gabaergic' 'glutamatergic' };
Glia = 'microglia astrocyte oligodendrocyte';
CortexN = 'pyramidal basket gabaergic glutamatergic';
CerebN  = 'purkinje granule'; 

Species = {'Mouse', 'Rat', 'Human', 'Monkey'};

%% Soma R
edges = 0:0.03:30;
GliaR = zeros(1,length(edges));
GliaSR = [];
ConR  = zeros(1,length(edges));
ConSR = [];
CenR  = zeros(1,length(edges));
CenSR = [];

% loop over cell types
for i = 1:length(total_stats)
    % loop over total_stats and identify species and cell type

    SomaR = total_stats(i).soma_Reff( ~isnan( total_stats(i).soma_Reff ) );
  
    if isempty(SomaR)
        continue
    end

    SomaR = SomaR(~isoutlier(SomaR));

    [S, ~] = ksdensity(SomaR, edges);
    % normalise
    S = S/sum(S);

    % select color for cell type
    if contains(Glia, total_stats(i).celltype)
        GliaR = GliaR + S;
        GliaSR = [GliaSR; SomaR];
    elseif contains(CortexN, total_stats(i).celltype)
        ConR = ConR + S;
        ConSR = [ConSR; SomaR];
    elseif contains(CerebN, total_stats(i).celltype)
        CenR = CenR + S;
        CenSR = [CenSR; SomaR];
    end

end

[GliaR, ~] = ksdensity(GliaSR, edges);
GliaR = GliaR/sum(GliaR);

[ConR, ~] = ksdensity(ConSR, edges);
ConR = ConR/sum(ConR);

[CenR, ~] = ksdensity(CenSR, edges);
CenR = CenR/sum(CenR);

figure
hold on
fill([edges,max(edges),0], [GliaR,0,0], [233, 60, 1]/255, 'FaceAlpha', 0.25);
fill([edges,max(edges),0], [ConR,0,0],  [0, 119, 194]/255, 'FaceAlpha', 0.25);
fill([edges,max(edges),0], [CenR,0,0],  [92, 169,1]/255, 'FaceAlpha', 0.25);

title('Soma Radius')
xlabel('R_{soma} (\mu)')
ylabel('pdf');
xlim([0 6])
%% Branch Length
edges = 1:0.1:300;
GliaL = zeros(1,length(edges));
ConL  = zeros(1,length(edges));
CenL  = zeros(1,length(edges));

for i = 1:length(total_stats)

    % loop over total_stats and identify species and cell type
    BranchL =  total_stats(i).mean_L_branch( ~isnan( total_stats(i).mean_L_branch ) );

    if isempty(BranchL) || std(BranchL)<1
        continue
    end

    BranchL = BranchL(~isoutlier(BranchL));

    [BL, x] = ksdensity(BranchL, edges);
    % normalise
    BL = BL/sum(BL);

    x = log10(x);
    ind = x>0;

    % select color for cell type
    if contains(Glia, total_stats(i).celltype)
        GliaL = GliaL + BL;
        colour = [233, 60, 1]/255;
    elseif contains(CortexN, total_stats(i).celltype)
        ConL = ConL + BL;
        colour = [0, 119, 194]/255;
    elseif contains(CerebN, total_stats(i).celltype)
        CenL = CenL + BL;
        colour = [92, 169,1]/255;
    end     

end
GliaL = GliaL/sum(GliaL);
ConL = ConL/sum(ConL);
CenL = CenL/sum(CenL);

figure
hold on
fill([log10(edges),max(log10(edges)),min(log10(edges))], [GliaL,0,0], [233, 60, 1]/255, 'FaceAlpha', 0.25);
fill([log10(edges),max(log10(edges)),min(log10(edges))], [ConL,0,0],  [0, 119, 194]/255, 'FaceAlpha', 0.25);
fill([log10(edges),max(log10(edges)),min(log10(edges))], [CenL,0,0],  [92, 169,1]/255, 'FaceAlpha', 0.25);

title('Branch length')
xticks([0 1 2 3])
xticklabels({'10^0','10^1','10^2','10^3'})
xlabel('Branch length (\mum)')
ylabel('pdf');
% xlim([0 20])


%% Intra cellular residence time
Permeability = 20 * 10^-3; % \mum-ms
edges = 0:0.25:200;

GliacellSV     = [];
GliacellSVmean = [];
GliaCellTex    = zeros(1,length(edges));

CortexNcellSV     = [];
CortexNcellSVmean = [];
CortexNCellTex    = zeros(1,length(edges));

CerebNcellSV     = [];
CerebNcellSVmean = [];
CerebNCellTex    = zeros(1,length(edges));
% loop over cell types
for i = 1:length(total_stats)
    % loop over total_stats and identify species and cell type

    S = total_stats(i).soma_surf + total_stats(i).branch_S;
    V = total_stats(i).soma_vol  + total_stats(i).branch_V;
    cellSV = S./V;
    
    if isempty(cellSV)
        continue
    end
    
    cellSV = cellSV(~isnan(cellSV));
    cellSV = cellSV(~isoutlier(cellSV));
    
    Tex = 1./(cellSV*Permeability);
    
    [Tex, ~] = ksdensity(Tex, edges);
    % normalise
    Tex = Tex/sum(Tex);  
    
    % select color for cell type
    if contains(Glia, total_stats(i).celltype)
        GliacellSV     = [GliacellSV; cellSV];
        GliacellSVmean = [GliacellSVmean; mean(cellSV)];
        GliaCellTex    = GliaCellTex + Tex;
    elseif contains(CortexN, total_stats(i).celltype)
        CortexNcellSV     = [CortexNcellSV; cellSV];
        CortexNcellSVmean = [CortexNcellSVmean; mean(cellSV)];
        CortexNCellTex    = CortexNCellTex + Tex;       
    elseif contains(CerebN, total_stats(i).celltype)
        CerebNcellSV     = [CerebNcellSV; cellSV];
        CerebNcellSVmean = [CerebNcellSVmean; mean(cellSV)];
        CerebNCellTex = CerebNCellTex + Tex;
    end
   
end

%% Intra branch residence time
Permeability = 20 * 10^-3; % \mum-ms
edges = 0:1:200;
edgesSV = 0:1:100;

GliabranchSV     = [];%zeros(1,length(edgesSV)-1);
GliabranchSVmean = [];
GliaBranchTex    = zeros(1,length(edges));
CortexNbranchSV     = [];%zeros(1,length(edgesSV)-1);
CortexNbranchSVmean = [];
CortexNBranchTex    = zeros(1,length(edges));
CerebNbranchSV     = [];%zeros(1,length(edgesSV)-1);
CerebNbranchSVmean = [];
CerebNBranchTex    = zeros(1,length(edges));

% loop over cell types
for i = 1:length(total_stats)

    cellSV = total_stats(i).mean_S_V_branch;
    
    if isempty(cellSV)
        continue
    end
      
    cellSV = cellSV(~isnan(cellSV));
    cellSV = cellSV(~isoutlier(cellSV));
    
    Tex = 1./(cellSV*Permeability);
    
    [Tex, ~] = ksdensity(Tex, edges);
    % normalise
    Tex = Tex/sum(Tex);  
    
    % select color for cell type
    if contains(Glia, total_stats(i).celltype)
        GliabranchSV     = [GliabranchSV; cellSV];
        GliabranchSVmean = [GliabranchSVmean; mean(cellSV)];
        GliaBranchTex    = GliaBranchTex + Tex;
        colour = [233, 60, 1]/255;
        
    elseif contains(CortexN, total_stats(i).celltype)
        CortexNbranchSV     = [CortexNbranchSV; cellSV];
        CortexNbranchSVmean = [CortexNbranchSVmean; mean(cellSV)];
        CortexNBranchTex    = CortexNBranchTex + Tex;
        colour = [0, 119, 194]/255;
        
    elseif contains(CerebN, total_stats(i).celltype)
        CerebNbranchSV     = [CerebNbranchSV; cellSV];
        CerebNbranchSVmean = [CerebNbranchSVmean; mean(cellSV)];
        CerebNBranchTex = CerebNBranchTex + Tex;
        colour = [92, 169, 1]/255;
    end
     
end

%% Cell residence & exchange time
figure 
subplot(1,3,1)

Permeability = 20*10^-3;
edges = 0:0.1:100;
TiG = 1./(GliacellSV*Permeability);
[TiG, x] = ksdensity(TiG, edges);
% normalise
TiG = TiG/sum(TiG); 

TiCo = 1./(CortexNcellSV*Permeability);
[TiCo, ~] = ksdensity(TiCo, edges);
% normalise
TiCo = TiCo/sum(TiCo); 

TiCe = 1./(CerebNcellSV*Permeability);
[TiCe, ~] = ksdensity(TiCe, edges);
% normalise
TiCe = TiCe/sum(TiCe);

subplot(1,4,1)
patch = fill([10*x,max(10*x),0], [TiG,0,0], [233, 60, 1 ]/255);
set(patch, 'FaceAlpha', 0.25);
hold on
patch = fill([10*x,max(10*x),0], [TiCo,0,0], [0, 119, 194]/255);
set(patch, 'FaceAlpha', 0.25);
patch = fill([10*x,max(10*x),0], [TiCe,0,0], [92, 169, 1 ]/255);
set(patch, 'FaceAlpha', 0.25);hold on
xlim([0 600])
title({'Intra cellular residence times', [ '\kappa = ' num2str(Permeability*100) '(\mumms^{-1})' ]})
xlabel('\fontsize{18}\tau_{i} \fontsize{12} (ms)')
ylabel('pdf')
axis square

subplot(1,4,2)
patch = fill([x,max(x),0], [TiG,0,0], [233, 60, 1 ]/255);
set(patch, 'FaceAlpha', 0.25);
hold on
patch = fill([x,max(x),0], [TiCo,0,0], [0, 119, 194]/255);
set(patch, 'FaceAlpha', 0.25);
patch = fill([x,max(x),0], [TiCe,0,0], [92, 169, 1 ]/255);
set(patch, 'FaceAlpha', 0.25);hold on
xlim([0 600])
title({'Intra cellular residence times', [ '\kappa = ' num2str(Permeability*1000) '(\mumms^{-1})' ]})
xlabel('\fontsize{18}\tau_{i} \fontsize{12} (ms)')
ylabel('pdf')
legend('Glia', 'Cortex neurons', 'Cerebelum neurons')
axis square


perm = (2:0.1:35)*10^-3;
fe = 0.3; % extracellular volume fraction

TiG  = 1./(GliacellSV*perm);    TiGM = mean(TiG,1);   TiGs = std(TiG,1);
TiCo = 1./(CortexNcellSV*perm); TiCoM = mean(TiCo,1); TiCos = std(TiCo,1);
TiCe = 1./(CerebNcellSV*perm);  TiCeM = mean(TiCe,1); TiCes = std(TiCe,1);

subplot(1,4,3)
plot(perm*1000, TiGM,  'Color', [233, 60, 1 ]/255, 'LineWidth', 1)
hold on
plot(perm*1000, TiCoM, 'Color', [0, 119, 194]/255, 'LineWidth', 1)
plot(perm*1000, TiCeM, 'Color', [92, 169, 1 ]/255, 'LineWidth', 1)

grid on
xlim([0 35])
ylim([0 max([TiGM, TiCoM, TiCeM])])
title('Intra cell residence time')
xlabel('Permeability (\mumms^{-1})')
ylabel('\fontsize{18}\tau_{i} \fontsize{12} (ms)')
axis square

TexG   = TiGM*fe; 
TexCo  = TiCoM*fe; 
TexCe  = TiCeM*fe;

subplot(1,4,4)
plot(perm*1000, TexG,  'Color', [233, 60, 1 ]/255, 'LineWidth', 1)
hold on
plot(perm*1000, TexCo, 'Color', [0, 119, 194]/255, 'LineWidth', 1)
plot(perm*1000, TexCe, 'Color', [92, 169, 1 ]/255, 'LineWidth', 1)

grid on
xlim([0 35])
ylim([0 max([TiGM, TiCoM, TiCeM])])
title('Cell exchange time')
xlabel('Permeability (\mumms^{-1})')
ylabel('\fontsize{18}\tau_{ex} \fontsize{12} (ms)')
legend('Glia', 'Cortex neurons', 'Cerebelum neurons')
axis square

%% Branch residence & exchange time
figure 
subplot(1,3,1)

Permeability = 20*10^-3;
edges = 0:0.1:100;
TiG = 1./(GliabranchSV*Permeability);
[TiG, x] = ksdensity(TiG, edges);
% normalise
TiG = TiG/sum(TiG); 

TiCo = 1./(CortexNbranchSV*Permeability);
[TiCo, ~] = ksdensity(TiCo, edges);
% normalise
TiCo = TiCo/sum(TiCo); 

TiCe = 1./(CerebNbranchSV*Permeability);
[TiCe, ~] = ksdensity(TiCe, edges);
% normalise
TiCe = TiCe/sum(TiCe);

subplot(1,4,1)
patch = fill([x*10,max(x*10),0], [TiG,0,0], [233, 60, 1 ]/255);
set(patch, 'FaceAlpha', 0.25);
hold on
patch = fill([x*10,max(x*10),0], [TiCo,0,0], [0, 119, 194]/255);
set(patch, 'FaceAlpha', 0.25);
patch = fill([x*10,max(x*10),0], [TiCe,0,0], [92, 169, 1 ]/255);
set(patch, 'FaceAlpha', 0.25);hold on
title({'Intra branch residence times', [ '\kappa = ' num2str(Permeability*100) '(\mumms^{-1})' ]})
xlabel('\fontsize{18}\tau_{i} \fontsize{12} (ms)')
ylabel('pdf')
axis square
xlim([0 600])

subplot(1,4,2)
patch = fill([x,max(x),0], [TiG,0,0], [233, 60, 1 ]/255);
set(patch, 'FaceAlpha', 0.25);
hold on
patch = fill([x,max(x),0], [TiCo,0,0], [0, 119, 194]/255);
set(patch, 'FaceAlpha', 0.25);
patch = fill([x,max(x),0], [TiCe,0,0], [92, 169, 1 ]/255);
set(patch, 'FaceAlpha', 0.25);hold on
title({'Intra branch residence times', [ '\kappa = ' num2str(Permeability*1000) '(\mumms^{-1})' ]})
xlabel('\fontsize{18}\tau_{i} \fontsize{12} (ms)')
ylabel('pdf')
legend('Glia', 'Cortex neurons', 'Cerebelum neurons')
axis square
xlim([0 600])


perm = (2:0.1:35)*10^-3;
fe = 0.3; % extracellular volume fraction

TiG  = 1./(GliabranchSVmean*perm);    TiGM = mean(TiG,1);   TiGs = std(TiG,1);
TiCo = 1./(CortexNbranchSVmean*perm); TiCoM = mean(TiCo,1); TiCos = std(TiCo,1);
TiCe = 1./(CerebNbranchSVmean*perm);  TiCeM = mean(TiCe,1); TiCes = std(TiCe,1);

subplot(1,4,3)
plot(perm*1000, TiGM,  'Color', [233, 60, 1 ]/255, 'LineWidth', 1)
hold on
plot(perm*1000, TiCoM, 'Color', [0, 119, 194]/255, 'LineWidth', 1)
plot(perm*1000, TiCeM, 'Color', [92, 169, 1 ]/255, 'LineWidth', 1)

grid on
xlim([0 35])
ylim([0 max([TiGM, TiCoM, TiCeM])])
title('Intra branch residence time')
xlabel('Permeability (\mumms^{-1})')
ylabel('\fontsize{18}\tau_{i} \fontsize{12}(ms)')
axis square


TexG   = TiGM*fe; 
TexCo  = TiCoM*fe; 
TexCe  = TiCeM*fe;

subplot(1,4,4)
plot(perm*1000, TexG,  'Color', [233, 60, 1 ]/255, 'LineWidth', 1)
hold on
plot(perm*1000, TexCo, 'Color', [0, 119, 194]/255, 'LineWidth', 1)
plot(perm*1000, TexCe, 'Color', [92, 169, 1 ]/255, 'LineWidth', 1)

grid on
xlim([0 35])
ylim([0 max([TiGM, TiCoM, TiCeM])])
title('Branch exchange time')
xlabel('Permeability (\mumms^{-1})')
ylabel('\fontsize{18}\tau_{ex} \fontsize{12}(ms)')
legend('Glia', 'Cortex neurons', 'Cerebelum neurons')
axis square
