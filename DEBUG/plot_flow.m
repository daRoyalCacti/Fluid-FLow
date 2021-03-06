%plots the flow for a give time
function [vx, vy, vz] = plot_flow(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f');
    fclose(fileID);
    x = data(1:5:end);
    y = data(2:5:end);
    
    vx = data(3:5:end);
    vy = data(4:5:end);
    vz = data(5:5:end);
    
%     figure
    quiver(x,y,vx,vy, 3)
    axis([min(x), max(x), min(y), max(y)])
%     quiver(y,x,vy,vx, 3)
end