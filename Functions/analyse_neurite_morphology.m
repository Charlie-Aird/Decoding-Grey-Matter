function [basic_stats, advanced_stats]=analyse_neurite_morphology(list_swc_files, outputfolder)


Ncells = numel(list_swc_files);

error = zeros(Ncells,1);
basic_stats = cell(Ncells,1);
advanced_stats = cell(Ncells,1);


for i = 1:Ncells
    
    cellpath = list_swc_files(i).folder;
    cellname = list_swc_files(i).name;
    
    tmp_dir_name = [cellpath '/tmp']; mkdir(tmp_dir_name);
    
    % Load SWC file
    try

    tmp = struct;
        
    tree = load_tree([cellpath '/' cellname]); 
    tree = repair_tree(tree);
    
    % Compute basic stats 
    basic_stats{i} = stats_tree(tree,[],[],'-x');
    
    % Compute advanced stats
    sect = dissect_tree(tree);
    
    bsect = sect(tree.R(sect(:,2))>1,:); % get starting and final node of a branch
    
    % Number of branches
    tmp.number_of_branches = size(bsect,1);
    
    % Surface, volume and surface to volume ration of each branch
    seg_surf = surf_tree(tree); % get segments surface
    seg_vol = vol_tree(tree); % get segments volume
    sL = len_tree(tree); % get segments length
    
    tmp.branch_surf = zeros(tmp.number_of_branches,1); 
    tmp.branch_vol = zeros(tmp.number_of_branches,1); 
    tmp.branch_surf_to_vol = zeros(tmp.number_of_branches,1);
    tmp.along_branch_D = cell(tmp.number_of_branches,1);
    tmp.branch_subsegments_L = cell(tmp.number_of_branches,1);
    tmp.branch_main_vec = zeros(tmp.number_of_branches,3);
    tmp.branch_subsegments_vec = cell(tmp.number_of_branches,1);
    tmp.branch_subsegments_angle = cell(tmp.number_of_branches,1);
    tmp.branch_linear_path = zeros(tmp.number_of_branches,6);
    tmp.branch_subsegments_path = cell(tmp.number_of_branches,1);
    tmp.branch_curve = cell(tmp.number_of_branches,1);
    tmp.branch_cumulative_curve_length = cell(tmp.number_of_branches,1);
    tmp.branch_curvature_radius = cell(tmp.number_of_branches,1);
    tmp.branch_curvature_vec = cell(tmp.number_of_branches,1);
    tmp.branch_microOD = zeros(tmp.number_of_branches,1);
    
    ipar = ipar_tree (tree);
    
    for ii=1:tmp.number_of_branches
        
        branch_node_idx = ipar(bsect (ii,2), ipar(bsect (ii,2),:)>=bsect (ii,1));
        
        tmp.branch_surf(ii) = sum(seg_surf( branch_node_idx ));
        tmp.branch_vol(ii) = sum(seg_vol( branch_node_idx ));
        tmp.branch_surf_to_vol(ii) = tmp.branch_surf(ii)./tmp.branch_vol(ii);
        tmp.along_branch_D{ii} = tree.D(branch_node_idx);
        tmp.branch_subsegments_L{ii} = sL(branch_node_idx);
        tmp.branch_main_vec(ii,:) = [tree.X(bsect(ii,2))-tree.X(bsect(ii,1)) tree.Y(bsect(ii,2))-tree.Y(bsect(ii,1)) tree.Z(bsect(ii,2))-tree.Z(bsect(ii,1))]; % get main branch vec
        tmp.branch_main_vec(ii,:) = tmp.branch_main_vec(ii,:)./norm(tmp.branch_main_vec(ii,:));
        tmp.branch_linear_path(ii,:) = [tree.X(bsect(ii,1)) tree.X(bsect(ii,2)) tree.Y(bsect(ii,1)) tree.Y(bsect(ii,2)) tree.Z(bsect(ii,1)) tree.Z(bsect(ii,2))];
        
        tmp.branch_curve{ii} = [];
        tmp.branch_subsegments_vec{ii} = [];
        tmp.branch_subsegments_angle{ii} = [];
        tmp.branch_subsegments_path{ii} = [];
        
        for kk = 1:numel(branch_node_idx)-1
            jj = branch_node_idx(kk);
            
            tmp.branch_curve{ii} = [tmp.branch_curve{ii}; tree.X(jj) tree.Y(jj) tree.Z(jj)];
            tmp_vec = [tree.X(jj+1)-tree.X(jj) tree.Y(jj+1)-tree.Y(jj) tree.Z(jj+1)-tree.Z(jj)];
            tmp_vec = tmp_vec./norm(tmp_vec);
            tmp.branch_subsegments_vec{ii} = [tmp.branch_subsegments_vec{ii}; tmp_vec]; % get subsegment vec
            tmp.branch_subsegments_angle{ii} = [tmp.branch_subsegments_angle{ii} acos(dot(tmp_vec,tmp.branch_main_vec(ii,:)))]; % get the angle between subsegment and main branch vec
            tmp.branch_subsegments_path{ii} = [tmp.branch_subsegments_path{ii}; tree.X(jj) tree.X(jj+1) tree.Y(jj) tree.Y(jj+1) tree.Z(jj) tree.Z(jj+1)];
        end
        tmp.branch_curve{ii} = [tmp.branch_curve{ii}; tree.X(bsect(ii,1)) tree.Y(bsect(ii,1)) tree.Z(bsect(ii,1))];
        [L,R,K] = curvature(tmp.branch_curve{ii});
        tmp.branch_cumulative_curve_length{ii} = L;
        tmp.branch_curvature_radius{ii} = R;
        tmp.branch_curvature_vec{ii} = K;
        tmp.branch_microOD(ii) = nanmean( sin(tmp.branch_subsegments_angle{ii}).^2 );
    end
    
    dist_from_soma = sqrt( (tree.X-tree.X(1)).^2 + (tree.Y-tree.Y(1)).^2 + (tree.Z-tree.Z(1)).^2);
    tmp.cell_domain_R = nanmean(max(dist_from_soma));
    
    advanced_stats{i} = tmp;
    
    catch
        error(i) = i;
    end
    
end

save([outputfolder '/neurites_stats.mat'], 'basic_stats', 'advanced_stats');

end