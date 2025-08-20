%% Meshing script

% location of cellular reconstructions
morphodata_path = '/cubric/data/c1441567/LargeScaleAnalysis2/data_for_TINS_review/morphological_analysis2/';
list = dir([morphodata_path '*']);
% output path
outPath = [pwd '/Meshes'];
if ~exist( outPath, 'dir' )
    mkdir( outPath );
end
% what species and cell type to mesh
species = "Human";%["Rat", "Mouse", "Monkey", "Human"];
celltypes = {'Microglia' 'Astrocyte' 'Oligodendrocyte' 'Pyramidal' 'Granule' 'Purkinje' 'Basket' 'Gabaergic' 'Glutamatergic' };
% number of meshes per cell type
n = 50;

% loop over species
for ff = 1:numel(species)
    % isolate folders for species
    list = dir([morphodata_path species{ff} '*']);
    % first level output folder
    outFolder1 = [outPath '/' species{ff}];
    if ~exist( outFolder1, 'dir' )
        mkdir(outFolder1)
    end
    % loop over corresponding driectories
    for i = numel(list)
        
        morphodata_filename = list(i).name;

        % find cell type
        for ii = 1:numel(celltypes)           
            if contains(morphodata_filename, celltypes{ii}, 'IgnoreCase',  true)
                celltype = celltypes{ii};
            end
        end
        % load trees
        load([morphodata_path list(i).name '/processed/trees.mat']);
        if isempty(trees)
            continue
        end
        % second level output folder
        outFolder2 = [outFolder1 '/' celltype];
        if ~exist( outFolder2, 'dir' )
            mkdir( outFolder2 );
        end
        % randomly select N trees
        nTrees = length(trees);
        N = min([nTrees, n]);       
        rN = randperm(nTrees, N);
        % loop over and mesh selected trees
        for ii = 1:N
            swc_tree(trees{rN(ii)}, [outFolder2 '/cell_' num2str(ii) '.swc'])
            blender_mesher_2(outFolder2, ['cell_' num2str(ii) '.swc'], max([0.05 min(trees{rN(ii)}.D/2)]), 1.4, 0.05 , ii);
        end
    end
end


