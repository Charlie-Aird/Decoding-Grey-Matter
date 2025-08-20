%% Analyse and clean cells

function [trees, error]=analyse_clean_cells(list_swc_files, outputfolder)

Ncells = numel(list_swc_files);

error = zeros(Ncells,1);

trees = cell(Ncells,1);
treesOG = cell(Ncells,1);

% loop over individual cells
for i = 1:Ncells
    % identify file path and name
    cellpath = list_swc_files(i).folder;
    cellname = list_swc_files(i).name;
    
    % skip if only representing the soma
    if contains(cellname, 'only_soma') || contains(cellname, 'full_cell')
        error(i) = 1;
        continue
    end
    
    % Load SWC file
    try
        tree = load_tree([cellpath '/' cellname]); 
        % some swc files contain multiple structures
        if iscell(tree)
            tree = tree{1};
        end
        % reject trees that are only composed of a few nodes
        if length(tree.R)<=10
            error(i) = 1;
            continue
        end
        tree = repair_tree(tree);
        
        % remove axonal component 
        comp = tree.R;
        uComp = unique( comp );  
        if contains(cellpath, "pyramidal") && length(uComp) == 4
            % pyramidal cells should be soma, basal, apical
            tree = delete_tree(tree, (comp==2) );
        elseif ~contains(cellpath, "pyramidal") && length(uComp) == 3
            % all other cells should be soma, densrite
            tree = delete_tree(tree, comp==2 );
        end 
        % remove if tree structure integrtiy comprimised or reduced to just the soma
        dA = tree.dA; dA = sum(dA,2);
        euc = eucl_tree(tree);
        insideSoma = euc <= tree.D(1)/2;
        if sum(dA==0)>=2 || sum(~insideSoma)==0
            tree = [];
        end
        trees{i} = tree;
    catch
        error(i) = 1;
        continue
    end
    
end
trees = trees(~cellfun('isempty', trees));
save([outputfolder '/trees.mat'], 'trees');

