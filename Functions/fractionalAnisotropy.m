
function [FAt, FAtS, Bingham, line] = fractionalAnisotropy(tree,l)

% reformat and resample the tree
tree = repair_tree(tree);
tree = reSoma_tree(tree);
tree = resample_tree(tree, l/2, 'none');
tree.D(1) = tree.D(2);

% deconstruct tree into constituent parts
[BLO, ~, ~]  = BLO_tree  (tree);
% fit segments to bars
count = 1;
% loop over bars
for i = 1:max(BLO)
%     try  %#ok<TRYNC>
    % create swc structure for each bar
    bar.X = tree.X(BLO==i);
    bar.Y = tree.Y(BLO==i);
    bar.Z = tree.Z(BLO==i);
    bar.D = tree.D(BLO==i);
    bar.R = tree.R(BLO==i);
    dA = tree.dA; dA(BLO~=i,:) = [];  dA(:,BLO~=i) = [];
    bar.dA = dA;
    bar.rnames = num2cell(1:max(bar.R));
    bar.name = num2str(i); 
    % decompose bars into sections of length l
    PL  = Pvec_tree(bar);
    PTS = [bar.X, bar.Y, bar.Z];
    euc = sqrt( sum( (PTS-PTS(1,:)).^2,2) );
    pl = ceil(euc/l);
    pl(pl==0) = 1;
    % loop over segments of length l
    while ~isempty(PTS) %j = 1:max(pl)
        % if the reamianing branch is too short move on
        if max(euc)<=0.75*l  
            break
        end
         
        ind = false(size(pl));
        for ii = 1:length(pl)
            if pl(ii)~=1
                break
            else
                ind(ii) = true;
            end
        end
        pts = PTS(ind,:);
        ptsM= mean(pts);
        
        if size(pts,1)==1
            pause(1)
        end
        
        % find line direction that fits points
        [~,~,v] = svd(pts-ptsM);
        % mean radius of points
        r = mean(bar.D(ind))/2;
        % save 
        line(count).pts = pts;
        line(count).u = v(:,1)';
        line(count).rSqr = r^2;
        line(count).comp = mode( bar.R(ind) );
        % increasec count number
        count = count + 1;       
        % reevaluate remaining points
        PL (ind) = [];
        PTS(ind,:) = [];
        if isempty(PTS)
            break
        end  
        euc = sqrt( sum( (PTS-PTS(1,:)).^2,2) );           
        pl = ceil(euc/l);
        pl(pl==0)=1;
    end
end
% compute the normalised line vectors
w = vertcat(line.rSqr);
w = sqrt(w/sum(w));
S = vertcat(line.u).*w;
% calculate orientation dispertion matrix T as scatter matrix of S
T = S'*S;
% find eigenvalues of T
tau = eig(T);
tauM= sum(tau)/3;
% if sum(tau)~=1
%     disp('eigenvalues not normalised')
% end
% compute fractional anisotropy of orientaion dispersion
FAt = sqrt( (3/2) * ( (tau(1) - tauM)^2 + (tau(2) - tauM)^2 + (tau(3) - tauM)^2 )/ sum(tau.^2));
% symetric fractional anisotropy
tauS = tau; tauS(1) = tauS(2); tauS = tauS/sum(tauS); tauSm = sum(tauS)/3;
FAtS = sqrt( (3/2) * ( (tauS(1) - tauSm)^2 + (tauS(2) - tauSm)^2 + (tauS(3) - tauSm)^2 )/ sum(tauS.^2));

%% Bingham/Watsom parameter estimation

vec = vertcat(line.u);
vec = [vec; -vec];

bb = BinghamDistribution.fit(vec', w);
Bingham = bb.Z;


    
    


