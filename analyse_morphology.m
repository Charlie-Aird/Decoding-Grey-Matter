%% Analyse Morphology 2

clear all
clc

diary('analyse_morphology_2_command_window_out.txt')

start_trees;
% addpath('/home/c1441567/Desktop/Code/Functions');
addpath('/cubric/data/c1441567/LargeScaleAnalysis/Functions');
blenderpath = '/Applications/Blender/blender.app/Contents/MacOS/blender';
addpath('/home/c1441567/MATLAB Add-Ons/Toolboxes/libDirectional-master')
startup

warning('off')

data_filename_list = {'Mouse_cells_for_TINS_review','Rat_cells_for_TINS_review', 'Monkey_cells_for_TINS_review', 'Human_cells_for_TINS_review'};
speciesType = ["Mouse"]; %["Mouse", "Rat", "Monkey", "Human"];
cellTypes = ["glutamatergic"];%["microglia", "astrocyte", "oligodendrocyte", "pyramidal", "granule", "basket", "purkinje", "gabaergic", "glutamatergic"];

% Length over which shape characteristics assessed
l = 10; % microns

numberOfCells = 0;

for ff = 1:numel(speciesType)
    
% find sub folders
raw_data_filename = data_filename_list{ff};
raw_data_path = '/cubric/data/c1441567/LargeScaleAnalysis2/data_for_TINS_review/';

raw_data_folder = [raw_data_path raw_data_filename];
% identify species
for j = 1:length(speciesType)
    if contains(raw_data_filename, speciesType{j})
        Species = speciesType{j}; 
    end 
end

% Read subdirectories within the main raw data directory
list = dir(raw_data_folder);

list_subdir=list(1);
k=1;
for i=1:numel(list)
    if list(i).isdir==1 && strcmp(list(i).name,'.')==0 && strcmp(list(i).name,'..')==0
        list_subdir(k) = list(i);
        k=k+1;
    end
end

% loop over file directories
for i = 1:numel(list_subdir) 
    
    % identify swc files
    folderpath = list_subdir(i).folder;
    foldername = list_subdir(i).name;
    list_swc_files = dir([folderpath '/' foldername '/**/*.swc']);
    % identify celltype
    for j = 1:length(cellTypes)
        if contains(foldername, cellTypes{j}, 'IgnoreCase', true)
            cellType = cellTypes{j};
            break
        else
            cellType = [];
            continue
        end 
    end 
    if isempty(cellType)
        continue
    end
  
    disp(['~~~~~~  ' Species ' - ' cellType '  ~~~~~~'])
    
    outputfolder = [raw_data_path 'morphological_analysis2/' raw_data_filename '_' foldername '/processed'];
    mkdir(outputfolder);
    
    % analyse and clean
    % look for saved trees
    if exist([outputfolder '/trees.mat'], 'file')==2
        disp('Loading trees')
        load([outputfolder '/trees.mat']);
        numberOfCells = numberOfCells + length(trees);
    else
    % otherwise load trees from swc files
    disp('Analysing trees')
    tic
    try 
        [trees, error] = analyse_clean_cells(list_swc_files, outputfolder, Species, cellType);
    catch
        disp('ERROR')
    end
    toc
    end
    
    % soma stats
    disp('Analysing soma')
    tic
    try
        cellpath = [raw_data_path '/Soma_Mesh' '/' Species '/' cellType];       
        soma_stats = analyse_soma_morphology_tree(trees, cellpath, outputfolder);
    catch
        disp('ERROR')
    end
    toc
    
    % neurite stats
    disp('Analysing neurite structure')
    tic
    try
        [basic_stats, advanced_stats] = analyse_neurite_morphology_tree(trees, outputfolder);
    catch
        disp('ERROR')
    end
    toc
    
    % topology
    disp('Analysing neurite topology')
    tic
    try
        [topology] = analyse_neurite_topology_tree(trees, outputfolder);
    catch
        disp('ERROR')
    end
    toc
    
    % shape
    disp('Analysing neurite shape')
    tic
    try
        shape = analyse_neurite_shape_tree(trees, l, outputfolder);
    catch
        disp('ERROR')
    end  
    toc
end
disp('Done')
end

diary('off')





