%% analyse Topology

function [topology_stat] =analyse_neurite_topology_tree(trees, outfolder)

Ncells = numel(trees);

topology_stat = cell(Ncells,1);
error = zeros(Ncells,1);

for i = 1:Ncells
    % identify and remove soma nodes (except origin)
    if length(unique(trees{i}.R))~=1
        delete = trees{i}.R==1;
        delete(1) = false;
        % delete soma nodes
        treeTMP = delete_tree(trees{i}, delete);
    else
        treeTMP = trees{i};
    end
    try
        % extract barcode
        barcode = barcode_tree(treeTMP);   
        topology_stat{i,1} = barcode;
    catch
        error(i) = 1;
    end
end
save([outfolder '/TopologyEuclideanDistance.mat'], 'topology_stat', 'error');