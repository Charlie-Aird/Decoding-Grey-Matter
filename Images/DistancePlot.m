%% bootsatrapping error
clear all

addpath '/cubric/data/c1441567/LargeScaleAnalysis2/functions'
morphodata_path = '/cubric/data/c1441567/LargeScaleAnalysis2/data_for_TINS_review/morphological_analysis2/';
list = dir([morphodata_path '*']);

%%
species = ["Rat", "Mouse", "Monkey", "Human"];
celltypes = {'Microglia' 'Astrocyte' 'Oligodendrocyte' 'Pyramidal' 'Granule' 'Purkinje' 'Basket' 'Gabaergic' 'Glutamatergic' };

BarCodes = struct;
edges = 0:100:1000;
[X, Y] = meshgrid(edges(1:end-1),edges(1:end-1));

% collate bar codes
% loop over species
for sp = 1:length(species)
    
    disp(['Species - ' species{sp}])
    
    % loop over cell types
    for ct = 1:length(celltypes)
        
        disp([' - ' celltypes{ct}])
        
        barcodes = [];
        BC = [];
        scaledbarcodes = cell(1,1);
        scaledBC = [];
        
        for i = 1:length(list)
            
            if contains(list(i).name, species{sp}) && contains(list(i).name, celltypes{ct})
                load([morphodata_path list(i).name '/processed/Topology.mat']);
                barcodes = [barcodes; topology_stat];
                BC = [BC; vertcat(topology_stat{:})]; 

                for ii = 1:length(topology_stat)
                    tt = abs( topology_stat{ii}./max(topology_stat{ii}(:,2)) );
                    tt(tt<1e-3) = 0;
                    scaledbarcodes{end+1,1} = tt;
                    scaledBC = [scaledBC; tt];
                end
            end            
        end
        
        if isempty(barcodes)
            continue
        end

        tic
        S = ksdensity(BC, [X(:) Y(:)]);
        pImage = reshape(S/sum(S(:)), size(X));
        toc
        
        BarCodes.( species{sp} ).( celltypes{ct} ).barcodes = barcodes;
        BarCodes.( species{sp} ).( celltypes{ct} ).BC       = BC;
        BarCodes.( species{sp} ).( celltypes{ct} ).pImage   = pImage;

        tic
        S = ksdensity(scaledBC, [sX(:) sY(:)]);
        SpImage = reshape(S/sum(S(:)), size(sX));
        toc
        
        BarCodes.( species{sp} ).( celltypes{ct} ).Sbarcodes = scaledbarcodes;
        BarCodes.( species{sp} ).( celltypes{ct} ).SBC       = scaledBC;
        BarCodes.( species{sp} ).( celltypes{ct} ).SpImage   = SpImage;
        
        % concatenate rat and mouse data
        if strcmp(species{sp}, 'Mouse')
            BarCodes.Rodent.( celltypes{ct} ).barcodes = [BarCodes.Rat.( celltypes{ct} ).barcodes; barcodes];
            BarCodes.Rodent.( celltypes{ct} ).BC       = [BarCodes.Rat.( celltypes{ct} ).BC; BC];
            S = ksdensity(BarCodes.( species{sp} ).( celltypes{ct} ).BC, [X(:) Y(:)]);
            pImage = reshape(S/sum(S(:)), size(X));
            BarCodes.Rodent.( celltypes{ct} ).pImage   = pImage;

            BarCodes.Rodent.( celltypes{ct} ).Sbarcodes = [BarCodes.Rat.( celltypes{ct} ).Sbarcodes; scaledbarcodes];
            BarCodes.Rodent.( celltypes{ct} ).SBC       = [BarCodes.Rat.( celltypes{ct} ).SBC; scaledBC];
            S = ksdensity(BarCodes.( species{sp} ).( celltypes{ct} ).SBC, [sX(:) sY(:)]);
            SpImage = reshape(S/sum(S(:)), size(sX));
            BarCodes.Rodent.( celltypes{ct} ).SpImage   = SpImage;
        end
        
        % concatenate rat and mouse data
        if strcmp(species{sp}, 'Human') || strcmp(species{sp}, 'Monkey')
            try
            BarCodes.Hominid.( celltypes{ct} ).barcodes = [BarCodes.Monkey.( celltypes{ct} ).barcodes; barcodes];
            BarCodes.Hominid.( celltypes{ct} ).BC       = [BarCodes.Monkey.( celltypes{ct} ).BC; BC];
            S = ksdensity(BarCodes.( species{sp} ).( celltypes{ct} ).BC, [X(:) Y(:)]);
            pImage = reshape(S/sum(S(:)), size(X));
            BarCodes.Hominid.( celltypes{ct} ).pImage   = pImage;
            catch
            BarCodes.Hominid.( celltypes{ct} ).barcodes = barcodes;
            BarCodes.Hominid.( celltypes{ct} ).BC       = BC;
            S = ksdensity(BarCodes.( species{sp} ).( celltypes{ct} ).BC, [sX(:) sY(:)]);
            pImage = reshape(S/sum(S(:)), size(sX));
            BarCodes.Hominid.( celltypes{ct} ).pImage   = pImage;
            end
        end
    end
end
%% Plot persitence images
figure
for i = 1:length(celltypes)
    subplot(3,3,i)
    imagesc(BarCodes.Rodent.(celltypes{i}).SpImage);
    axis image
    xticks([1 length(BarCodes.Rodent.(celltypes{i}).SpImage)])
    xticklabels({'0', '1'})
    yticks([1 length(BarCodes.Rodent.(celltypes{i}).SpImage)])
    yticklabels({'0', '1'})
    title(celltypes{i});
end

%% Plot distance Matrix
animal = 'Human';
celltypes = fields( BarCodes.( animal ) );

D  = zeros(length(celltypes));
SD = zeros(length(celltypes));
for i = 1:length(celltypes)
    
    A = BarCodes.( animal ).( celltypes{i} ).SpImage;
    indA = A~=0;
    
    for ii = i+1:length(celltypes)
        
        disp([celltypes{i} ' - ' celltypes{ii}])

        B = BarCodes.( animal ).( celltypes{ii} ).SpImage;
        indB = B~=0;

        D(i,ii) = sum( abs(A(:)-B(:)) );
        D(ii,i) = D(i,ii);
            
    end
end

figure
imagesc(D); axis image
xticks(1:length(celltypes))
xticklabels(celltypes)
yticks(1:length(celltypes))
yticklabels(celltypes)
colorbar
colormap turbo
caxis([0 2])

%% Compare between species

Species1 = 'Rodent';
cellTypes1 = fields(BarCodes.( Species1 ));
Species2 = 'Hominid';
cellTypes2 = fields(BarCodes.( Species2 ));
cellTypes = cellTypes1(contains(cellTypes1, cellTypes2)); 
cellType = ["Microglia", "Astrocyte", "Oligodendrocyte", "Pyramidal", "Granule", "Purkinje", "Basket", "Gabaergic", "Glutamatergic"];
cellType = cellType(contains(cellType, cellTypes));
DistanceM  = zeros(length(cellType));

for i = 1:length(cellType)
    A = BarCodes.(Species1).(cellType{i}).pImage;
    
    for ii = 1:length(cellType)
        B = BarCodes.(Species2).(cellType{ii}).pImage;
        D = (A - B);
        DistanceM (i,ii) = sum( abs(D(:)) );
    end
   
end
figure
imagesc(DistanceM); shading flat; axis image; colormap turbo; colorbar
xticks(1:length(D))
xticklabels(cellType)
yticks(1:length(D))
yticklabels(cellType)
colorbar
colormap turbo
caxis([0 2])

