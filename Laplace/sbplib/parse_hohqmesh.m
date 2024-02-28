function [Nnodes,Nedges,Nelements,Bound_order,coords,edge_data,element_data,bnames] = parse_hohqmesh(mesh_filename)

fileID = fopen(mesh_filename,'r');
A = fscanf(fileID,"%s\n",[1,1]);

B = fscanf(fileID,"          %d          %d          %d          %d\n",[4,1]);

Nnodes = B(1);
Nedges = B(2);
Nelements = B(3);
Bound_order = B(4);

all_bname = {};

for idx = 1:Nnodes
    coords(idx,:) = fscanf(fileID," %f      %f        %f     ",[1,3]);
end

%remove z
coords =  coords(:,1:2);

for idx = 1:Nedges
    edge_data(idx,:) = fscanf(fileID,"     %d       %d       %d       %d       %d       %d\n",[1,6]);
end

for idx = 1:Nelements
    element_data_tmp = fscanf(fileID,"           %d           %d           %d           %d\n",[1,4]);
    is_phys_boundary = fscanf(fileID,"           %d           %d           %d           %d\n",[1,4]);
    GL_points = [];
    for i = 1:4
        if is_phys_boundary(i) == 0
            continue
        end
        for GL_point_idx = 1:Bound_order+1
            GL_points_tmp = fscanf(fileID," %f       %f   %f     \n",[1,3]);


            % remove z
            GL_points_tmp =  GL_points_tmp(:,1:2);
            GL_points = [GL_points;GL_points_tmp];
        end
    end
    bound_name1 = fscanf(fileID," %s",1);
    all_bname{end+1} = char(bound_name1);

    bound_name2 = fscanf(fileID," %s",1);
    all_bname{end+1} = char(bound_name2);

    bound_name3 = fscanf(fileID," %s",1);
    all_bname{end+1} = char(bound_name3);

    bound_name4 = fscanf(fileID," %s",1);
    all_bname{end+1} = char(bound_name4);
    
    element_data{idx} = {element_data_tmp,is_phys_boundary,GL_points,{bound_name1,bound_name2,bound_name3,bound_name4}};
    fscanf(fileID," --- --- --- --- \n");
end

bnames = unique(all_bname,'stable');
bnames = bnames(2:end);