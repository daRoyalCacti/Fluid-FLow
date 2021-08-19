%plots the point cloud for a give time
function Plot_mesh(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f');
    fclose(fileID);
    x = data(1:3:end);
    y = data(2:3:end);
    z = data(3:3:end);
    
    figure
    plot3(x,y,z, 'o')
end