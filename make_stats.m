clear all
clc

%% Read all morpho data

celltypes = {'microglia' 'astrocyte' 'oligodendrocyte' 'pyramidal' 'granule' 'purkinje' 'basket' 'gabaergic' 'glutamatergic' };

morphodata_filename_list = {'Mouse','Rat', 'Monkey', 'Human'};
morphodata_path = '/cubric/data/c1441567/LargeScaleAnalysis2/data_for_TINS_review/morphological_analysis2/';

output_folder = 'complete_summary_stats';
mkdir(output_folder)

%% Create total stat

total_stats = struct;
total_shape = struct;
count = 1;
for ff = 1:numel(morphodata_filename_list)
    
    list = dir([morphodata_path morphodata_filename_list{ff} '*']);
    
    % For each animal store all the morphodata belonging to the same
    % celltype in the same variable
    
    for i = 1:numel(list)
        
        total_stats(count).animal = morphodata_filename_list{ff};
        total_shape(count).animal = morphodata_filename_list{ff};
        
        morphodata_filename = list(i).name;
        
        % find cell type
        for ii = 1:numel(celltypes)           
            if contains(morphodata_filename, celltypes{ii}, 'IgnoreCase',  true)
                total_stats(count).celltype = celltypes{ii};
                total_shape(count).celltype = celltypes{ii};
            end
        end
        
        file_input = [morphodata_path morphodata_filename '/processed'];
        % look in folder for soma data
        if exist([file_input '/soma_stats.mat'],'file')
            morphodata_soma = load([file_input '/soma_stats.mat'], 'soma_stats');
            soma_stats = morphodata_soma.soma_stats;
        else
            soma_stats = cell(1,1);
        end
        % look in folder for neurite data
        if exist([file_input '/neurites_stats.mat'],'file')
            morphodata_neurite = load([file_input '/neurites_stats.mat'], 'basic_stats', 'advanced_stats');
            basic_neurites_stats = morphodata_neurite.basic_stats;
            advanced_neurites_stats = morphodata_neurite.advanced_stats;
        else
            basic_neurites_stats = cell(1,1);
            advanced_neurites_stats = cell(1,1);
        end
        % look in folder for shape data
        if exist([file_input '/Shape.mat'],'file')
            morphodata_shape = load([file_input '/Shape.mat'], 'shape_stat');
            shape_stats = morphodata_shape.shape_stat;           
        else
            shape_stats = cell(1,1);
        end
        
        
        n_cells_soma = numel(soma_stats);
        
        total_stats(count).n_cells = n_cells_soma;
        total_stats(count).soma_vol = zeros(n_cells_soma, 1);
        total_stats(count).soma_surf = zeros(n_cells_soma, 1);
        total_stats(count).soma_Reff = zeros(n_cells_soma, 1);
        total_stats(count).soma_nominalR = zeros(n_cells_soma, 1);
        total_stats(count).soma_connectivityratio = zeros(n_cells_soma, 1);
        total_stats(count).n_projections = zeros(n_cells_soma, 1);
        
        n_cells_neurite = numel(advanced_neurites_stats);
        
        total_stats(count).mean_BO = zeros(n_cells_neurite, 1);
        total_stats(count).std_BO = zeros(n_cells_neurite, 1);
        total_stats(count).celldomain_R = zeros(n_cells_neurite, 1);
        total_stats(count).branch_S = zeros(n_cells_neurite, 1);
        total_stats(count).branch_V = zeros(n_cells_neurite, 1);
        total_stats(count).mean_S_V_branch = zeros(n_cells_neurite, 1);
        total_stats(count).mean_L_branch = zeros(n_cells_neurite, 1);
        total_stats(count).std_L_branch = zeros(n_cells_neurite, 1);
        total_stats(count).mean_mean_R_along_branch = zeros(n_cells_neurite, 1);
        total_stats(count).mean_std_R_along_branch = zeros(n_cells_neurite, 1);
        total_stats(count).coeff_V_R_along_branch = zeros(n_cells_neurite, 1);
        total_stats(count).mean_microOD_branch = zeros(n_cells_neurite, 1);
        total_stats(count).mean_Rcurvature_branch = zeros(n_cells_neurite, 1);
        total_stats(count).mean_tortuosity_branch = zeros(n_cells_neurite, 1);
        
        for ii = 1:n_cells_soma
            
            % Soma morphometrics
            
            if ~isempty(soma_stats{ii})
                
            total_stats(count).soma_vol(ii) = soma_stats{ii}.vol; 
            total_stats(count).soma_surf(ii) = soma_stats{ii}.surf;
            total_stats(count).soma_Reff(ii) = (3./(4.*pi).*soma_stats{ii}.vol).^(1/3);
            total_stats(count).soma_nominalR(ii) = soma_stats{ii}.nominal_D./2;
            total_stats(count).soma_connectivityratio(ii) = soma_stats{ii}.projsurf_over_somasurf; if total_stats(count).soma_connectivityratio(ii)>1, total_stats(count).soma_connectivityratio(ii)=1; end
            total_stats(count).n_projections(ii) = soma_stats{ii}.n_proj;
            
            end
            
        end
        
        for ii = 1:n_cells_neurite
            % Neurites morphometrics
            
            if ~isempty(basic_neurites_stats{ii}) && ~isempty(advanced_neurites_stats{ii})
                
                
                total_stats(count).mean_BO(ii) = nanmean(basic_neurites_stats{ii}.dstats.BO{1,1});
                total_stats(count).std_BO(ii) = nanstd(basic_neurites_stats{ii}.dstats.BO{1,1});
                total_stats(count).celldomain_R(ii) = advanced_neurites_stats{ii}.cell_domain_R;
                total_stats(count).branch_S(ii) = sum(advanced_neurites_stats{ii}.branch_surf);
                total_stats(count).branch_V(ii) = sum(advanced_neurites_stats{ii}.branch_vol);
                total_stats(count).mean_S_V_branch(ii) = nanmean(advanced_neurites_stats{ii}.branch_surf_to_vol);
                
                L_branch = [];
                mean_R_along_branch = [];
                std_R_along_branch = [];
                R_curve = [];
                empt = false(numel(advanced_neurites_stats{ii}.along_branch_D),1);
                for bb = 1:numel(advanced_neurites_stats{ii}.along_branch_D)
                    try
                        L_branch = [L_branch; advanced_neurites_stats{ii}.branch_cumulative_curve_length{bb}(end)];
                        mean_R_along_branch = [mean_R_along_branch; nanmean(advanced_neurites_stats{ii}.along_branch_D{bb}/2)];
                        std_R_along_branch = [std_R_along_branch; nanstd(advanced_neurites_stats{ii}.along_branch_D{bb})];
                        R_curve = [R_curve; advanced_neurites_stats{ii}.branch_curvature_radius{bb}]; R_curve(R_curve>5000) = nan; R_curve(R_curve==0) = nan;
                        
                    catch
                        empt(bb) = true;
                    end
                end
                
                total_stats(count).mean_L_branch(ii) = nanmean(L_branch);
                total_stats(count).std_L_branch(ii) = nanstd(L_branch);
                total_stats(count).mean_mean_R_along_branch(ii) = nanmean(mean_R_along_branch);
                total_stats(count).mean_std_R_along_branch(ii) = nanmean(std_R_along_branch);
                total_stats(count).coeff_V_R_along_branch(ii)  = nanmean(std_R_along_branch./mean_R_along_branch)*100;                
                total_stats(count).mean_microOD_branch(ii) = nanmean(advanced_neurites_stats{ii}.branch_microOD);
                total_stats(count).mean_Rcurvature_branch(ii) = nanmean(R_curve);
                
                direct_L_branch =  sqrt( (advanced_neurites_stats{ii}.branch_linear_path(~empt,2)-advanced_neurites_stats{ii}.branch_linear_path(~empt,1)).^2 + (advanced_neurites_stats{ii}.branch_linear_path(~empt,4)-advanced_neurites_stats{ii}.branch_linear_path(~empt,3)).^2 + (advanced_neurites_stats{ii}.branch_linear_path(~empt,6)-advanced_neurites_stats{ii}.branch_linear_path(~empt,5)).^2 );
                total_stats(count).mean_tortuosity_branch(ii) = nanmean(direct_L_branch./L_branch);
            end

        end
        
        count = count+1;
        
    end
     
end

% save([output_folder '/total_stats.mat'], 'total_stats');

%% Merge duplicate celltypes and generate final stats

% collate glial stats
stats_glia = struct;
stats_glia.soma_SV   = [];
stats_glia.soma_Reff = [];
stats_glia.soma_RmrD = [];
stats_glia.soma_connectivityratio = [];
stats_glia.n_projections = [];
stats_glia.mean_BO = [];
stats_glia.SVdomain = [];
stats_glia.celldomain_R = [];
stats_glia.mean_S_V_branch = [];
stats_glia.mean_L_branch = [];
stats_glia.mean_R_along_branch = [];
stats_glia.Reffective = [];
stats_glia.coeff_V_R_along_branch = [];
stats_glia.mean_microOD_branch = [];
stats_glia.mean_Rcurvature_branch = [];
stats_glia.mean_tortuosity_branch = [];
% collate neuron stats
stats_neuron= struct;
stats_neuron.soma_SV   = [];
stats_neuron.soma_Reff = [];
stats_neuron.soma_RmrD = [];
stats_neuron.soma_connectivityratio = [];
stats_neuron.n_projections = [];
stats_neuron.mean_BO = [];
stats_neuron.SVdomain = [];
stats_neuron.celldomain_R = [];
stats_neuron.mean_S_V_branch = [];
stats_neuron.mean_L_branch = [];
stats_neuron.mean_R_along_branch = [];
stats_neuron.Reffective = [];
stats_neuron.coeff_V_R_along_branch = [];
stats_neuron.mean_microOD_branch = [];
stats_neuron.mean_Rcurvature_branch = [];
stats_neuron.mean_tortuosity_branch = [];
% 
Y = [];
YY = [];
for j = 1:numel(celltypes) % loop over the celltypes
    
    celltype = celltypes{j};
         
    for  ff = 1:numel(morphodata_filename_list) % loop over the animals
        
        animal = morphodata_filename_list{ff};
        
        % group rat/mouse
        if contains('MouseRat', animal)
            animal = 'MouseRat';
        end    
        
        stats = struct;
        
        stats.animal = animal;
        stats.celltype = celltype;
        
        stats.soma_vol  = [];
        stats.soma_surf = [];
        stats.soma_SV   = [];
        stats.soma_Reff = [];
        stats.soma_nominalR = [];
        stats.soma_connectivityratio = [];
        
        stats.n_projections = [];
        stats.mean_BO = [];
        stats.std_BO = [];
        stats.SVdomain = [];
        stats.celldomain_R = [];
        stats.mean_S_V_branch = [];
        stats.mean_L_branch = [];
        stats.std_L_branch = [];
        stats.mean_mean_R_along_branch = [];
        stats.mean_std_R_along_branch = [];
        stats.coeff_V_R_along_branch = [];
        stats.mean_microOD_branch = [];
        stats.mean_Rcurvature_branch = [];
        stats.mean_tortuosity_branch = [];

        stats.glia = struct;
        
        for i = 1:numel(total_stats)
            % group stats
            if contains(animal, total_stats(i).animal) && strcmp(celltype, total_stats(i).celltype)
                
                % General stats
                stats.celldomain_R  = [stats.celldomain_R; total_stats(i).celldomain_R]; stats.celldomain_R(stats.celldomain_R==0) = nan;
                stats.n_projections = [stats.n_projections; total_stats(i).n_projections]; stats.n_projections(stats.n_projections==0) = nan;
                stats.mean_BO       = [stats.mean_BO; total_stats(i).mean_BO]; stats.mean_BO(stats.mean_BO==0) = nan;
                stats.std_BO        = [stats.std_BO; total_stats(i).std_BO]; stats.std_BO(stats.std_BO==0) = nan;
                stats.SVdomain      = [stats.SVdomain; (total_stats(i).soma_surf + total_stats(i).branch_S)./(total_stats(i).soma_vol + total_stats(i).branch_V)]; stats.SVdomain(stats.SVdomain==0) = nan;
                
                % Soma stats
                stats.soma_vol  = [stats.soma_vol; total_stats(i).soma_vol]; stats.soma_vol(stats.soma_vol<1e-5) = nan;
                stats.soma_surf = [stats.soma_surf; total_stats(i).soma_surf]; stats.soma_surf(stats.soma_surf==0) = nan;
                stats.soma_SV   = [stats.soma_SV; total_stats(i).soma_surf./total_stats(i).soma_vol]; stats.soma_SV(stats.soma_SV==0) = nan; stats.soma_SV(stats.soma_SV>1e+2) = nan;
                stats.soma_Reff = [stats.soma_Reff; total_stats(i).soma_Reff]; stats.soma_Reff(stats.soma_Reff==0) = nan;
                stats.soma_nominalR = [stats.soma_nominalR; total_stats(i).soma_nominalR]; stats.soma_nominalR(stats.soma_nominalR==0) = nan;
                stats.soma_connectivityratio = [stats.soma_connectivityratio; total_stats(i).soma_connectivityratio]; stats.soma_connectivityratio(stats.soma_connectivityratio==0) = nan;                
                
                % Branch stats
                stats.mean_S_V_branch = [stats.mean_S_V_branch; total_stats(i).mean_S_V_branch]; stats.mean_S_V_branch(stats.mean_S_V_branch==0) = nan; stats.soma_SV(stats.mean_S_V_branch==inf) = nan;
                stats.mean_L_branch   = [stats.mean_L_branch; total_stats(i).mean_L_branch]; stats.mean_L_branch(stats.mean_L_branch==0) = nan;
                stats.std_L_branch    = [stats.std_L_branch; total_stats(i).std_L_branch]; stats.std_L_branch(stats.std_L_branch==0) = nan;
                stats.mean_mean_R_along_branch = [stats.mean_mean_R_along_branch; total_stats(i).mean_mean_R_along_branch]; stats.mean_mean_R_along_branch(stats.mean_mean_R_along_branch==0) = nan;
                stats.mean_std_R_along_branch  = [stats.mean_std_R_along_branch; total_stats(i).mean_std_R_along_branch]; stats.mean_std_R_along_branch(stats.mean_std_R_along_branch==0) = nan;
                stats.coeff_V_R_along_branch   = [stats.coeff_V_R_along_branch; total_stats(i).coeff_V_R_along_branch]; stats.coeff_V_R_along_branch(stats.coeff_V_R_along_branch==0) = nan;
                stats.mean_microOD_branch      = [stats.mean_microOD_branch; total_stats(i).mean_microOD_branch]; stats.mean_microOD_branch(stats.mean_microOD_branch==0) = nan;
                stats.mean_Rcurvature_branch   = [stats.mean_Rcurvature_branch; total_stats(i).mean_Rcurvature_branch]; stats.mean_Rcurvature_branch(stats.mean_Rcurvature_branch==0) = nan; stats.mean_Rcurvature_branch(stats.mean_Rcurvature_branch==inf) = nan;
                stats.mean_tortuosity_branch   = [stats.mean_tortuosity_branch; total_stats(i).mean_tortuosity_branch]; stats.mean_tortuosity_branch(stats.mean_tortuosity_branch==0) = nan;

            end  
            
        end
                    
        stats_summary = struct;

        stats_summary.animal = animal;
        stats_summary.celltype = celltype;
        stats_summary.metrics = {'mean' 'standard_deviation' '1st_quartile', 'median', '3rd_quartile'};
        
        stats_summary.Ncells_with_soma_meas = numel(stats.soma_vol);
        stats_summary.Ncells_with_neurites_meas = numel(stats.mean_L_branch);
        % General
        stats_summary.celldomain_R  = [nanmean(stats.celldomain_R)  nanstd(stats.celldomain_R) quantile(stats.celldomain_R,[0.25 0.50 0.75])];
        stats_summary.n_projections = [nanmean(stats.n_projections) nanstd(stats.n_projections) quantile(stats.n_projections,[0.25 0.50 0.75])];
        stats_summary.mean_BO       = [nanmean(stats.mean_BO)       nanstd(stats.mean_BO) quantile(stats.mean_BO,[0.25 0.50 0.75])];
        stats_summary.SVdomain      = [nanmean(stats.SVdomain)      nanstd(stats.SVdomain) quantile(stats.SVdomain,[0.25 0.50 0.75])];
        
        % Soma
%         stats_summary.soma_vol  = [nanmean(stats.soma_vol) nanstd(stats.soma_vol) quantile(stats.soma_vol,[0.25 0.50 0.75])];
%         stats_summary.soma_surf = [nanmean(stats.soma_surf) nanstd(stats.soma_surf) quantile(stats.soma_surf,[0.25 0.50 0.75])];        
        stats_summary.soma_SV   = [nanmean(stats.soma_surf./stats.soma_vol) nanstd(stats.soma_surf./stats.soma_vol) quantile(stats.soma_surf./stats.soma_vol,[0.25 0.50 0.75])];      
        stats_summary.soma_Reff = [nanmean(stats.soma_Reff) nanstd(stats.soma_Reff) quantile(stats.soma_Reff,[0.25 0.50 0.75])];        
        stats_summary.soma_RmrD = ( nanmean( stats.soma_Reff( ~isoutlier(stats.soma_Reff) ).^7 ) / nanmean( stats.soma_Reff( ~isoutlier(stats.soma_Reff) ).^3 ) ).^(1/4) ;
        stats_summary.soma_nominalR = [nanmean(stats.soma_nominalR) nanstd(stats.soma_nominalR) quantile(stats.soma_nominalR,[0.25 0.50 0.75])];
        stats_summary.soma_connectivityratio = [nanmean(stats.soma_connectivityratio) nanstd(stats.soma_connectivityratio) quantile(stats.soma_connectivityratio,[0.25 0.50 0.75])];

        % Branch
        stats_summary.mean_S_V_branch = [nanmean(stats.mean_S_V_branch) nanstd(stats.mean_S_V_branch) quantile(stats.mean_S_V_branch,[0.25 0.50 0.75])];
        stats_summary.mean_L_branch   = [nanmean(stats.mean_L_branch) nanstd(stats.mean_L_branch) quantile(stats.mean_L_branch,[0.25 0.50 0.75])];
        stats_summary.mean_mean_R_along_branch = [nanmean(stats.mean_mean_R_along_branch) nanstd(stats.mean_mean_R_along_branch) quantile(stats.mean_mean_R_along_branch,[0.25 0.50 0.75])];
        stats_summary.Reffective = [ ( nanmean(stats.mean_mean_R_along_branch( ~isoutlier(stats.mean_mean_R_along_branch) ).^6 ) / nanmean(stats.mean_mean_R_along_branch( ~isoutlier(stats.mean_mean_R_along_branch) ).^2) ).^(1/4) ];
        stats_summary.coeff_V_R_along_branch   = [nanmean(stats.coeff_V_R_along_branch) nanstd(stats.coeff_V_R_along_branch) quantile(stats.coeff_V_R_along_branch,[0.25 0.50 0.75])];
        stats_summary.mean_microOD_branch = [nanmean(stats.mean_microOD_branch) nanstd(stats.mean_microOD_branch) quantile(stats.mean_microOD_branch,[0.25 0.50 0.75])];
        stats_summary.mean_Rcurvature_branch = [nanmean(stats.mean_Rcurvature_branch) nanstd(stats.mean_Rcurvature_branch) quantile(stats.mean_Rcurvature_branch,[0.25 0.50 0.75])];
        stats_summary.mean_tortuosity_branch = [nanmean(stats.mean_tortuosity_branch) nanstd(stats.mean_tortuosity_branch) quantile(stats.mean_tortuosity_branch,[0.25 0.50 0.75])];

        save([output_folder '/total_soma_stats_' animal '_' celltype '.mat'], 'stats', 'stats_summary')
        
        strct = stats_summary;
        txt = fieldnames(strct) ;
        X = struct2cell(strct); 
        Y = [Y; X(1:length(X)~=3)'];
        X = cell2table(X);       
        T = [txt X];
        writetable(T, [output_folder '/summary_stats_' animal '_' celltype '1.xlsx']);

        % glia
        if contains( 'microglia_astrocyte_oligodendrocyte', stats.celltype)
            try
            stats_glia.soma_SV   = [stats_glia.soma_SV; [min(stats.soma_SV(stats.soma_SV<100)), max(stats.soma_SV(stats.soma_SV<100)), nanmean(stats.soma_SV(stats.soma_SV<100))]];
            stats_glia.soma_Reff = [stats_glia.soma_Reff; [min(stats.soma_Reff(~isoutlier(stats.soma_Reff))), max(stats.soma_Reff(~isoutlier(stats.soma_Reff))), nanmean(stats.soma_Reff(~isoutlier(stats.soma_Reff)))]];
            stats_glia.soma_RmrD = [stats_glia.soma_RmrD; [stats_summary.soma_RmrD, stats_summary.soma_RmrD, stats_summary.soma_RmrD]];
            stats_glia.soma_connectivityratio = [stats_glia.soma_connectivityratio; [min(stats.soma_connectivityratio(~isoutlier(stats.soma_connectivityratio))), max(stats.soma_connectivityratio(~isoutlier(stats.soma_connectivityratio))), nanmean(stats.soma_connectivityratio(~isoutlier(stats.soma_connectivityratio)))]];
            
            stats_glia.n_projections = [stats_glia.n_projections; [min(stats.n_projections), max(stats.n_projections), nanmean(stats.n_projections)]];
            stats_glia.mean_BO = [stats_glia.mean_BO; [min(stats.mean_BO), max(stats.mean_BO), nanmean(stats.mean_BO)]];
            stats_glia.SVdomain = [stats_glia.SVdomain; [min(stats.SVdomain), max(stats.SVdomain), nanmean(stats.SVdomain)]];
            stats_glia.celldomain_R = [stats_glia.celldomain_R; [min(stats.celldomain_R(~isoutlier(stats.celldomain_R))), max(stats.celldomain_R(~isoutlier(stats.celldomain_R))), nanmean(stats.celldomain_R(~isoutlier(stats.celldomain_R)))]];
            
            stats_glia.mean_S_V_branch = [stats_glia.mean_S_V_branch; [min(stats.mean_S_V_branch), max(stats.mean_S_V_branch), nanmean(stats.mean_S_V_branch)]];
            stats_glia.mean_L_branch = [stats_glia.mean_L_branch; [min(stats.mean_L_branch(~isoutlier(stats.mean_L_branch))), max(stats.mean_L_branch(~isoutlier(stats.mean_L_branch))), nanmean(stats.mean_L_branch(~isoutlier(stats.mean_L_branch)))]];
            stats_glia.mean_R_along_branch = [stats_glia.mean_R_along_branch; [min(stats.mean_mean_R_along_branch(~isoutlier(stats.mean_mean_R_along_branch))), max(stats.mean_mean_R_along_branch(~isoutlier(stats.mean_mean_R_along_branch))), nanmean(stats.mean_mean_R_along_branch(~isoutlier(stats.mean_mean_R_along_branch)))]];
            stats_glia.Reffective = [stats_glia.Reffective; [ stats_summary.Reffective, stats_summary.Reffective, stats_summary.Reffective]];
            stats_glia.coeff_V_R_along_branch = [stats_glia.coeff_V_R_along_branch; [min(stats.coeff_V_R_along_branch(~isoutlier(stats.coeff_V_R_along_branch))), max(stats.coeff_V_R_along_branch(~isoutlier(stats.coeff_V_R_along_branch))), nanmean(stats.coeff_V_R_along_branch(~isoutlier(stats.coeff_V_R_along_branch)))]];
            stats_glia.mean_microOD_branch = [stats_glia.mean_microOD_branch; [min(stats.mean_microOD_branch), max(stats.mean_microOD_branch), nanmean(stats.mean_microOD_branch)]];
            stats_glia.mean_Rcurvature_branch = [stats_glia.mean_Rcurvature_branch; [min(stats.mean_Rcurvature_branch), max(stats.mean_Rcurvature_branch), nanmean(stats.mean_Rcurvature_branch)]];
            stats_glia.mean_tortuosity_branch = [stats_glia.mean_tortuosity_branch; [min(stats.mean_tortuosity_branch), max(stats.mean_tortuosity_branch), nanmean(stats.mean_tortuosity_branch)]];
            end
        else % neuronal cells
            try
            stats_neuron.soma_SV   = [stats_neuron.soma_SV; [min(stats.soma_SV), max(stats.soma_SV), nanmean(stats.soma_SV)]];
            stats_neuron.soma_Reff = [stats_neuron.soma_Reff; [min(stats.soma_Reff(~isoutlier(stats.soma_Reff))), max(stats.soma_Reff(~isoutlier(stats.soma_Reff))), nanmean(stats.soma_Reff(~isoutlier(stats.soma_Reff)))]];
            stats_neuron.soma_RmrD = [stats_neuron.soma_RmrD; [stats_summary.soma_RmrD, stats_summary.soma_RmrD, stats_summary.soma_RmrD]];
            stats_neuron.soma_connectivityratio = [stats_neuron.soma_connectivityratio; [min(stats.soma_connectivityratio(~isoutlier(stats.soma_connectivityratio))), max(stats.soma_connectivityratio(~isoutlier(stats.soma_connectivityratio))), nanmean(stats.soma_connectivityratio(~isoutlier(stats.soma_connectivityratio)))]];
            
            stats_neuron.n_projections = [stats_neuron.n_projections; [min(stats.n_projections), max(stats.n_projections), nanmean(stats.n_projections)]];
            stats_neuron.mean_BO = [stats_neuron.mean_BO; [min(stats.mean_BO), max(stats.mean_BO), nanmean(stats.mean_BO)]];
            stats_neuron.SVdomain = [stats_neuron.SVdomain; [min(stats.SVdomain), max(stats.SVdomain), nanmean(stats.SVdomain)]];
            stats_neuron.celldomain_R = [stats_neuron.celldomain_R; [min(stats.celldomain_R(~isoutlier(stats.celldomain_R))), max(stats.celldomain_R(~isoutlier(stats.celldomain_R))), nanmean(stats.celldomain_R(~isoutlier(stats.celldomain_R)))]];
            
            stats_neuron.mean_S_V_branch = [stats_neuron.mean_S_V_branch; [min(stats.mean_S_V_branch), max(stats.mean_S_V_branch), nanmean(stats.mean_S_V_branch)]];
            stats_neuron.mean_L_branch = [stats_neuron.mean_L_branch; [min(stats.mean_L_branch(~isoutlier(stats.mean_L_branch))), max(stats.mean_L_branch(~isoutlier(stats.mean_L_branch))), nanmean(stats.mean_L_branch(~isoutlier(stats.mean_L_branch)))]];
            stats_neuron.mean_R_along_branch = [stats_neuron.mean_R_along_branch; [min(stats.mean_mean_R_along_branch(~isoutlier(stats.mean_mean_R_along_branch))), max(stats.mean_mean_R_along_branch(~isoutlier(stats.mean_mean_R_along_branch))), nanmean(stats.mean_mean_R_along_branch(~isoutlier(stats.mean_mean_R_along_branch)))]];
            stats_neuron.Reffective = [stats_neuron.Reffective; [ stats_summary.Reffective, stats_summary.Reffective, stats_summary.Reffective]];
            stats_neuron.coeff_V_R_along_branch = [stats_neuron.coeff_V_R_along_branch; [min(stats.coeff_V_R_along_branch), max(stats.coeff_V_R_along_branch), nanmean(stats.coeff_V_R_along_branch)]];
            stats_neuron.mean_microOD_branch = [stats_neuron.mean_microOD_branch; [min(stats.mean_microOD_branch), max(stats.mean_microOD_branch), nanmean(stats.mean_microOD_branch)]];
            stats_neuron.mean_Rcurvature_branch = [stats_neuron.mean_Rcurvature_branch; [min(stats.mean_Rcurvature_branch), max(stats.mean_Rcurvature_branch), nanmean(stats.mean_Rcurvature_branch)]];
            stats_neuron.mean_tortuosity_branch = [stats_neuron.mean_tortuosity_branch; [min(stats.mean_tortuosity_branch), max(stats.mean_tortuosity_branch), nanmean(stats.mean_tortuosity_branch)]];
            end
        end
        
    end

end

%%
metrics = repmat(X.X(3), 1, size(Y,2));
ind = cellfun(@length, Y);
ind = ind~=5;
ind = ind(1,:);
ind(1) = 1;
metrics(ind) = {'Mean'};

TT = [txt(1:length(txt)~=3)'; metrics; Y];
TT = cell2table(TT);

writetable(TT, [output_folder '/summary_stats_full.xlsx']);

sTT = size(TT);
for col = 1:sTT(1)
    for row = 4:sTT(2)
        if isnumeric( TT{col,row}{1} )
            TT{col,row}{1} = round( TT{col,row}{1}, 2, 'significant' );
        end
    end
end
               
writetable(TT, [output_folder '/summary_stats_full_rounded.xlsx']);

%
metrics = repmat(XX.XX(3), 1, size(YY,2));
ind = cellfun(@length, YY);
ind = ind~=5;
ind = ind(1,:);
ind(1) = 1;
metrics(ind) = {'Mean'};

TS = [txts(1:length(txts)~=3)'; metrics; YY];
TS = cell2table(TS);

writetable(TS, [output_folder '/summary_shape_full.xlsx']);

sTS = size(TS);
for col = 1:sTS(1)
    for row = 4:sTS(2)
        if isnumeric( TS{col,row}{1} )
            TS{col,row}{1} = round( TS{col,row}{1}, 2, 'significant' );
        end
    end
end
               
writetable(TS, [output_folder '/summary_shape_full_rounded.xlsx']);