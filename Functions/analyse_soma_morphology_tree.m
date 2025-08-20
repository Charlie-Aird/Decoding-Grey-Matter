function [soma_stats]=analyse_soma_morphology_tree(trees, cellpath, outputfolder)

Ncells = numel(trees);

error = zeros(Ncells,1);
soma_stats = cell(Ncells,1);

for i = 1:Ncells

    tmp_dir_name = [cellpath '/tmp']; mkdir(tmp_dir_name);
    
    % Load SWC file
    try
        
    tree = trees{i};
    
    soma_nominal_D = nanmean(tree.D(tree.R==1));
    first_proj = dist_tree(tree,soma_nominal_D/2*1.1);

    n_proj = sum(full(first_proj)==1);
    D_proj = tree.D(full(first_proj)==1);
    dist_trh = soma_nominal_D./2.*1.0;
    
    % Identify soma and closest projections
    
    %rtree = resample_tree(tree, 5);
    rtree = tree;
    within_soma = sqrt((rtree.X-rtree.X(1)).^2 + (rtree.Y-rtree.Y(1)).^2 + (rtree.Z-rtree.Z(1)).^2)<=dist_trh;

    only_soma_tree = rtree;

    only_soma_tree.X(within_soma==0) = [];
    only_soma_tree.Y(within_soma==0) = [];
    only_soma_tree.Z(within_soma==0) = [];tree = trees{i};
    
    soma_nominal_D = nanmean(tree.D(tree.R==1));
    first_proj = dist_tree(tree,soma_nominal_D/2*1.1);

    n_proj = sum(full(first_proj)==1);
    D_proj = tree.D(full(first_proj)==1);
    dist_trh = soma_nominal_D./2.*1.0;
    
    % Identify soma and closest projections
    
    %rtree = resample_tree(tree, 5);
    rtree = tree;
    within_soma = sqrt((rtree.X-rtree.X(1)).^2 + (rtree.Y-rtree.Y(1)).^2 + (rtree.Z-rtree.Z(1)).^2)<=dist_trh;

    only_soma_tree = rtree;

    only_soma_tree.X(within_soma==0) = [];
    only_soma_tree.Y(within_soma==0) = [];
    only_soma_tree.Z(within_soma==0) = [];
    only_soma_tree.D(within_soma==0) = [];
    only_soma_tree.R(within_soma==0) = [];
    
    idx = find(within_soma==1);
    
    dA = full(only_soma_tree.dA);
    
    dAnew = dA(idx,idx); 
    only_soma_tree.dA = sparse(dAnew);
    
    % Save SWC of only soma

    swc_tree(only_soma_tree,[tmp_dir_name '/only_soma.swc']);

    % Mesh only soma and compute soma volume and surface area
    
    [vertex, face, vol, surf] = blender_mesher_stats(tmp_dir_name, 'only_soma.swc', 0.1, i);
    
%     h = figure; hold on
%     options.curvature_smoothing = 10;
%     options.verb = 0;
%     [~,~,~,~,~,Cgauss,~] = compute_curvature(vertex,face,options);
%     %options.face_vertex_color = perform_saturation(Cgauss,1.2);
%     options.face_vertex_color = Cgauss;
%     plot_mesh(vertex,face, options); shading interp; colormap jet(256); colorbar, axis([-20 20 -20 20 -20 20]), hold on, plot3(linspace(0,10,10)', ones(10,1).*-20, ones(10,1).*(0), 'k-','LineWidth', 3.0),  view(2)
%     title('Soma shape and Gaussian curvature')
%     saveas(h,[outputfolder '/soma_' cellname '.png']);
%     close(h)

    soma_stats{i}.nominal_D = soma_nominal_D;
    soma_stats{i}.n_proj = n_proj;
    soma_stats{i}.D_proj = D_proj;

    soma_stats{i}.vertex = vertex;
    soma_stats{i}.face = face;
    soma_stats{i}.vol = vol;
    soma_stats{i}.surf = surf;
    soma_stats{i}.error = error(i);
    
    soma_stats{i}.projsurf_over_somasurf = sum(pi.*(D_proj./2).^2)./surf;
    
    catch
        error(i) = i;
    end  
end

save([outputfolder '/soma_stats.mat'], 'soma_stats', 'error');