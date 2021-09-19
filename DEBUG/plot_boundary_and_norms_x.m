    file_loc = "boundary_points.txt";
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
    file_loc = "normal_vectors.txt";
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f');
    fclose(fileID);
    y = data(1:5:end);
    z = data(2:5:end);
    
    vx = data(3:5:end);
    vy = data(4:5:end);
    vz = data(5:5:end);
    
    hold on
    quiver(y,z,vy,vz)
    
    has_x = vx~=0;
    plot(y(has_x), z(has_x), '.')
    