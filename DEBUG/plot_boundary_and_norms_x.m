    file_loc = "boundary_points_x.txt";
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f');
    fclose(fileID);
    y = data(1:3:end);
    z = data(2:3:end);
    
    b = data(3:3:end);
    
    ys = y(b==1);
    zs = z(b==1);
    figure
    plot(ys, zs, 'o')
    
%plots the flow for a give time
    file_loc = "normal_vectors_x.txt";
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f');
    fclose(fileID);
    y = data(1:5:end);
    z = data(2:5:end);
    
    vx = data(3:5:end);
    vy = data(4:5:end);
    vz = data(5:5:end);
    
    hold on
    div = (max(z) - min(z) ) / (max(y) - min(y));
    quiver(y,z,vy/div,vz, 'ShowArrowHead', 'off')
    
    has_x = vx~=0;
    plot(y(has_x), z(has_x), '.')
    daspect([1 1 1])
% daspect([1 40 1])
    xlabel('y')
    ylabel('z')
    