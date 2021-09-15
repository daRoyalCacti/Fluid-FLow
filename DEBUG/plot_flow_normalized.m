%plots the flow for a give time
function [vx, vy, vz] = plot_flow_normalized(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f');
    fclose(fileID);
    x = data(1:5:end);
    y = data(2:5:end);
    
    vx = data(3:5:end);
    vy = data(4:5:end);
    vz = data(5:5:end);
    
%     figure
    quiver(x,y, vx./sqrt(vx.^2+vy.^2), vy./sqrt(vx.^2+vy.^2))
    axis([min(x), max(x), min(y), max(y)])
end