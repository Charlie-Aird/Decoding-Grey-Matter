%% analyse shape

function shape_stat = analyse_neurite_shape_tree(trees, l, outfolder)

Ncells = numel(trees);

shape = struct();
shape_stat = cell(Ncells,1);
error =false(Ncells,1);

for i = 1:Ncells  
    disp(['Cell - '  num2str(i) ])
%     try       
    % fractional anistropy
    [shape.FAt, shape.FAtS, shape.Bingham, shape.line] = fractionalAnisotropy(trees{i}, l);
    
    shape_stat{i,1} = shape;
%     catch
%         error(i) = 1;
%     end
end
save([outfolder '/Shape.mat'], "shape_stat", "error");
 