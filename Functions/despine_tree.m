% despine_tree    Cumulative summation along termial branches of a tree.
% 
% ds_intree = despine_tree (intree, v, options)
% -------------------------------------
%
% Identifies Termial points and their associated branches
% Saves the nodes that the branch is composed of as well as total segment length 
% If branch length is bellow a length and volume threshold then the branch
% is removed
% Threshold can be defined by user but if not then length threshold is 8 (um)
% relating to the larger spine length as stated in 
% 'Sampling issues in quantitative analysis of dendritic spines morphology'
%
% Uses       len_tree; vol_tree; delete_tree

function ds_intree = despine_tree (intree,max_spine_length,max_spine_volume)

% trees : contains the tree structures in the trees package
global       trees

if (nargin < 1) || isempty (intree)
    % {DEFAULT tree: last tree in trees cell array}
    intree   = length (trees);
end

if nargin >= 1 && nargin < 3
    max_spine_length = 3;
    max_spine_volume = 10;
end

% Adgaceny Matrix
dA = intree.dA;
dA = full(dA);

% Identify Terminal Points
T = find((ones (1, size (dA, 1)) * dA) == 0)';
% Identify Segment lengths and volumes
V = vol_tree(intree);
L = len_tree(intree);

% Record Terminal Branch Nodes in 'Branch'
% Record Terminal Node Length in Branch_length
Branch = cell(length(T),1);
Branch_length = zeros(1,length(T));
Branch_parent = zeros(1,length(T));
Branch_volume = zeros(1,length(T));
  
for i = 1:length(T)
    
%     Branch{i} = [];
    
    % Identify Connected Nodes
    c = T(i);    
    % Level of connectedness (is zero for terminal nodes)
    con = 0;
    
    % Iterate Until Segment Reaches Branching Point (con=2)
    while con < 2
        
        Branch{i} = [Branch{i}; c];
        Branch_length(i) = Branch_length(i) + L(c);
        Branch_volume(i) = Branch_volume(i) + V(c);
        
        % Identify Connected Nodes
        c = find(dA(c,:));
        % Connections on Previous Node
        con = sum(dA(:,c)); 
    end
    Branch_parent(i) = c;
end

% find terminal points that share an immediate parent
[C,~,ic] = unique(Branch_parent);
shared_p = C(histcounts(ic,1:numel(C)+1) >1 );

ignore = [];
% compare branch lengths of branches with mutual parent 
for p = shared_p

    % identify terminal nodes with mutual parent
    t = T(Branch_parent==p);
    % idnetify the longer branch
    ind = Branch_length(Branch_parent==p)>min(Branch_length(Branch_parent==p));
    % ignore the longer branch
    ignore = [ignore; find(T==t(ind))];   
    
end

Branch(ignore) = [];
Branch_length(ignore) = [];
Branch_volume(ignore) = [];

% Find index (idx) of spines given they are smaller than mean Branch_length and some
% volume
idx = Branch_length<max_spine_length & Branch_volume<max_spine_volume ;
ds_intree = delete_tree(intree, cell2mat(Branch(idx)));

