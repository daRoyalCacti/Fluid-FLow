    file_loc = "boundary_points.txt";
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f');
    fclose(fileID);
    x = data(1:3:end);
    y = data(2:3:end);
    
    b = data(3:3:end);
    
    xs = x(b==1);
    ys = y(b==1);
    figure
    plot(xs, ys, 'o')
    
%plots the flow for a give time
    file_loc = "normal_vectors.txt";
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f');
    fclose(fileID);
    x = data(1:5:end);
    y = data(2:5:end);
    
    vx = data(3:5:end);
    vy = data(4:5:end);
    vz = data(5:5:end);
    
    hold on
    quiver(x,y,vx,vy)
    
    has_z = vz~=0;
    plot(x(has_z), y(has_z), '.')