%plots the point cloud for a give time
function Plot_mesh(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f %f %f %f %f');
    fclose(fileID);
    x = data(1:9:end);
    y = data(2:9:end);
    z = data(3:9:end);
    
    vx = data(4:9:end);
    vy = data(5:9:end);
    vz = data(6:9:end);
    
    nx = data(7:9:end);
    ny = data(8:9:end);
    nz = data(9:9:end);
    
%     figure
    plot3(x,y,z, 'o')
    %hold on
    %quiver3(x,y,z, vx, vy, vz, 'r')
    %quiver3(x,y,z, nx, ny, nz, 'm')
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
end