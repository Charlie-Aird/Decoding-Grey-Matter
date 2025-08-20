% BARCODE_TREE   Barcodes from persistence homology.
% (trees package)
%
% barcode = barcode_tree (intree, options)
% ----------------------------------------
%
% Returns the 
%
% Input
% -----
% - intree   ::integer:  index of tree in trees or structured tree
% - options  ::string:
%     '-2d'  : 2-dimensional lengths / Not implemented yet
%     '-s'   : show / Not implemented yet
%     {DEFAULT: ''}
%
% Output
% -------
% - barcode  ::Nx1 vector: length values of each segment
%
% Example
% -------
% len_tree     (sample_tree, '-s')
%
% 
% Interpreted from Kanari et al 2018 with help from Lida Kanari.
%
% See also BLO_tree, 
% Uses BLO_tree
%
% the TREES toolbox: edit, generate, visualise and analyse neuronal trees
% Copyright (C) 2009 - 2023  Hermann Cuntz

function [barcode, diameter] = barcode_tree (intree, V, options)

% trees : contains the tree structures in the trees package
global       trees

ver_tree     (intree); % verify that input is a tree structure

if (nargin < 1) || isempty (intree)
    % {DEFAULT tree: last tree in trees cell array}
    intree   = length (trees);
end

ver_tree     (intree); % verify that input is a tree structure

% use full tree for this function
if ~isstruct (intree)
    tree     = trees {intree};
else
    tree     = intree;
end

if (nargin < 2) || isempty (V)
    % {DEFAULT: length of segment in um, this way path length determines
    % order}
    V        = len_tree (tree);
end

if (nargin < 3) || isempty (options)
    % {DEFAULT: no option}
    options  = '';
end

[BLO, ~, clBLO]  = BLO_tree  (tree, V);
pV               = Pvec_tree (tree, V);
N                = max (BLO);
D                = tree.D;
birth            = zeros (N, 1);
death            = zeros (N, 1);
diameter = cell(N,1);
for counter      = 1 : N
    death (counter) = max (pV (BLO == counter));
    birth (counter) = death (counter) - max (clBLO (BLO == counter));
    
    diameter{counter} = resample( D(BLO == counter), 100, sum(BLO == counter) );
end

barcode          = [birth death];



