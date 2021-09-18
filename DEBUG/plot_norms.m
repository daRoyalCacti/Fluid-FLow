%plots the flow for a give time
% function plot_norms()
    file_loc = "normal_vectors.txt";
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f');
    fclose(fileID);
    x = data(1:5:end);
    y = data(2:5:end);
    
    vx = data(3:5:end);
    vy = data(4:5:end);
    %vz = data(5:5:end);
    
    figure
    quiver(x,y,vx,vy)
% end