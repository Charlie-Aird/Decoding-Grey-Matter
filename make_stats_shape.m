clear all
clc

%% Read all morpho data

celltypes = {'microglia' 'astrocyte' 'oligodendrocyte' 'pyramidal' 'granule' 'purkinje' 'basket' 'gabaergic' 'glutamatergic' };

morphodata_filename_list = {'Mouse','Rat', 'Monkey', 'Human'};
morphodata_path = '/cubric/data/c1441567/LargeScaleAnalysis2/data_for_TINS_review/morphological_analysis2/';

output_folder = 'complete_summary_stats';
mkdir(output_folder)

%% Create total stat

total_shape = struct;
count = 1;
for ff = 1:numel(morphodata_filename_list)
    
    list = dir([morphodata_path morphodata_filename_list{ff} '*']);
    
    % For each animal store all the morphodata belonging to the same
    % celltype in the same variable
    
    for i = 1:numel(list)
        
        total_shape(count).animal = morphodata_filename_list{ff};
        
        morphodata_filename = list(i).name;
        
        % find cell type
        for ii = 1:numel(celltypes)           
            if contains(morphodata_filename, celltypes{ii}, 'IgnoreCase',  true)
                total_shape(count).celltype = celltypes{ii};
            end
        end
        
        file_input = [morphodata_path morphodata_filename '/processed'];

        % look in folder for shape data
        if exist([file_input '/Shape.mat'],'file')
            morphodata_shape = load([file_input '/Shape.mat'], 'shape_stat');
            shape_stats = morphodata_shape.shape_stat;           
        else
            shape_stats = cell(1,1);
        end
        
        % shape stats
        n_cells_shape = length(shape_stats);
        total_shape(count).Bing  = zeros(n_cells_shape, 2);
        total_shape(count).Watson = zeros(n_cells_shape, 1); % Watson parameter
        total_shape(count).OD    = zeros(n_cells_shape, 1); % Orientational deispersion
        total_shape(count).tau   = zeros(n_cells_shape, 3); % eigenvalues
        total_shape(count).FA    = zeros(n_cells_shape, 1); % Fractional anisotropy
        total_shape(count).adjFA = zeros(n_cells_shape, 1); % adjusted FA

        for ii = 1:n_cells_shape
            
            line = shape_stats{ii}.line;
            
            vec = vertcat(line.u);
            vec = [vec; -vec];
            
            w = vertcat(line.rSqr);
            w = sqrt(w/sum(w));          
            w = [w; w];
            w = w/sum(w);

            bb = BinghamDistribution.fit(vec', w');
            Bingham = bb.Z;
            % find if cell is 2D
            if abs( Bingham(1) ) > abs( 4*Bingham(2) )
                kappa = abs( Bingham(2) );
                if kappa>10
                    continue
                end
            else
                kappa = abs( mean( Bingham(1:2) ) );
                if kappa>10
                    continue
                end
            end
            total_shape(count).Bing(ii,:)  = abs(Bingham(1:2));
            total_shape(count).Watson(ii) = kappa;
            total_shape(count).OD(ii)     = (2/pi)*atan( 1/kappa );
            
            % compute the normalised line vectors
            w = vertcat(shape_stats{ii, 1}.line.rSqr);
            w = sqrt(w/sum(w));
            S = vertcat(shape_stats{ii, 1}.line.u).*w;
            % calculate orientation dispertion matrix T as scatter matrix of S
            T = S'*S;
            % find eigenvalues of T
            tau = eig(T);
            tau = tau/sum(tau);
            
            total_shape(count).tau(ii,:) = tau;
            total_shape(count).FA(ii)    = shape_stats{ii}.FAt;
            total_shape(count).adjFA(ii) = shape_stats{ii}.FAtS;

        end
        
        count = count+1;
        
    end
     
end

%% Merge duplicate celltypes and generate final stats
 
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

        shape_summary = struct;
        shape_summary.animal = animal;
        shape_summary.celltype = celltype;
        shape_summary.metrics = {'mean' 'standard_deviation' '1st_quartile', 'median', '3rd_quartile'};
        
        stats_shape = struct;
        stats_shape.Watson = [];
        stats_shape.OD     = [];
        stats_shape.tau    = [];
        stats_shape.FA     = [];
        stats_shape.adjFA  = [];
        for i = 1:length(total_shape)
            if contains(animal, total_shape(i).animal) && strcmp(celltype, total_shape(i).celltype)
                stats_shape.Watson = [stats_shape.Watson; total_shape(i).Watson]; stats_shape.Watson(stats_shape.Watson==0) = nan;
                stats_shape.OD     = [stats_shape.OD; total_shape(i).OD]; stats_shape.Watson(stats_shape.OD==0) = nan;
                stats_shape.tau    = [stats_shape.tau; total_shape(i).tau]; stats_shape.tau(stats_shape.tau==0) = nan;
                stats_shape.FA     = [stats_shape.FA; total_shape(i).FA]; stats_shape.FA(stats_shape.FA==0) = nan;
                stats_shape.adjFA  = [stats_shape.adjFA; total_shape(i).adjFA]; stats_shape.adjFA(stats_shape.adjFA==0) = nan;

            end
        end
        
        shape_summary.Watson = [nanmean(stats_shape.Watson) nanstd(stats_shape.Watson) quantile(stats_shape.Watson, [0.25 0.50 0.75])];
        shape_summary.OD     = [nanmean(stats_shape.OD) nanstd(stats_shape.OD) quantile(stats_shape.OD, [0.25 0.50 0.75])];
        shape_summary.tau    = [nanmean(stats_shape.tau,1) nanstd(stats_shape.tau,1) ];%quantile(stats_shape.tau, [0.25 0.50 0.75], 1)];
        shape_summary.FA     = [nanmean(stats_shape.FA) nanstd(stats_shape.FA) quantile(stats_shape.FA, [0.25 0.50 0.75])];
        shape_summary.adjFA  = [nanmean(stats_shape.adjFA) nanstd(stats_shape.adjFA) quantile(stats_shape.adjFA, [0.25 0.50 0.75])];
        
        strct = shape_summary;
        txts = fieldnames(strct) ;
        XX = struct2cell(strct); 
        YY = [YY; XX(1:length(XX)~=3)'];
        XX = cell2table(XX);       
        TS = [txts XX];
        
    end

end

%%

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