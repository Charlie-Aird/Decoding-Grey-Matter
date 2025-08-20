function [stats] = blender_mesher_2(foldername, swc_filename, resolution, scaleFactor, decimate, i, Stats)

 %% Make the mesh (wrapper for Blender functions)+

[~,name,~] = fileparts(swc_filename); 
 
fid = fopen([foldername '/Bmesher' num2str(i) '.py'],'w');
 
fprintf(fid,'# How to use: blender -b -P mesher.py\n');
fprintf(fid,'import argparse\n');
fprintf(fid,'import bpy\n');
fprintf(fid,'import bmesh\n');
fprintf(fid,'filename = "%s.swc"\n',[foldername '/' name]); 
fprintf(fid,'bpy.data.scenes["Scene"].make_neuron_meta.neuron_file_name=filename\n');
fprintf(fid,'bpy.ops.mnm.make_line_mesh()\n');
fprintf(fid,'smallestRad = bpy.data.scenes["Scene"].make_neuron_meta.smallest_radius_in_file\n');
fprintf(fid,'largestRad = bpy.data.scenes["Scene"].make_neuron_meta.largest_radius_in_file\n');
%fprintf(fid,'meshres = min(0.6*smallestRad, 0.2*largestRad)\n');
%fprintf(fid,'meshres = 0.5*smallestRad\n');
fprintf(fid,['meshres =' num2str(resolution) '\n']);
fprintf(fid,'bpy.data.scenes["Scene"].make_neuron_meta.mesh_resolution=meshres\n');
fprintf(fid,['bpy.data.scenes["Scene"].make_neuron_meta.meta_ball_scale_factor=' num2str(scaleFactor) '\n']); % need to expand as meshing shrinks the diamter for some reason
%fprintf(fid,'bpy.data.scenes["Scene"].make_neuron_meta.meta_ball_threshold=0.0001\n');
fprintf(fid,'bpy.ops.mnm.make_neuron_from_data()\n');
%fprintf(fid,'bpy.data.metaballs["neuron"].threshold=0.01\n');
fprintf(fid,'bpy.ops.object.convert(target=''MESH'')\n');
fprintf(fid,'bpy.ops.object.join()\n');
fprintf(fid,'bpy.ops.object.modifier_add(type=''DECIMATE'')\n');
fprintf(fid,'namepostfix = "_cable_model"\n');
fprintf(fid,'SkeletonToBuild = "%s" + namepostfix\n', name);
fprintf(fid,'bpy.data.objects[SkeletonToBuild].modifiers["Decimate"].use_collapse_triangulate = 1\n');
fprintf(fid,['bpy.data.objects[SkeletonToBuild].modifiers["Decimate"].ratio =' num2str(decimate) '\n']); % triangle desimation reduce number of triangles
fprintf(fid,'bpy.ops.object.modifier_apply(apply_as=''DATA'')\n');
fprintf(fid,'bpy.ops.object.mode_set(mode=''EDIT'')\n');
fprintf(fid,'bpy.ops.object.modifier_add(type=''TRIANGULATE'')\n'); 
fprintf(fid,'bpy.ops.mesh.select_all(action=''TOGGLE'')\n');
fprintf(fid,'bpy.ops.mesh.quads_convert_to_tris(quad_method=''BEAUTY'')\n');
% stat calculations
fprintf(fid,'bpy.ops.mesh.print3d_info_volume()\n');
fprintf(fid,'bpy.ops.mesh.print3d_info_area()\n');
fprintf(fid,'bm = bmesh.new()\n');
fprintf(fid,'bm.from_object(bpy.context.object, bpy.context.scene)\n');
% volume
fprintf(fid,'print("VOLUMEi")\n');
fprintf(fid,'print(bm.calc_volume())\n'); 
fprintf(fid,'print("VOLUMEf")\n');
% surface area
fprintf(fid,'print("SURFACEAi")\n');
fprintf(fid,'print(sum(f.calc_area() for f in bm.faces))\n');
fprintf(fid,'print("SURFACEAf")\n');
% export
fprintf(fid,'bpy.ops.export_mesh.ply(filepath = "%s.ply")\n',[foldername '/' name]);
fprintf(fid,'bpy.ops.object.delete()\n');
 
fclose(fid);

command = ['/home/c1441567/Desktop/blender-2.78c/blender -b -P ' foldername '/Bmesher' num2str(i) '.py'];
[~,y] = system(command);
writing = strfind(y,'writing TEST');
disp( y(writing+12:end) )

% retrive stats
if (nargin > 6) && Stats
    vi = strfind(y,'VOLUMEi');
    vf = strfind(y,'VOLUMEf');
    ai = strfind(y,'SURFACEAi');
    af = strfind(y,'SURFACEAf');

    stats.vol = str2num(y(vi+7:vf-2));
    stats.surf = str2num(y(ai+9:af-2));
    
    save([swc_filename(1:end-4) '_stats.mat'], 'stats');     
end

delete([foldername '/Bmesher' num2str(i) '.py'])

end
