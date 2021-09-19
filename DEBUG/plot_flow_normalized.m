%plots the flow for a give time
function [vx, vy, vz,x,y] = plot_flow_normalized(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f');
    fclose(fileID);
    x = data(1:5:end);
    y = data(2:5:end);
    
    vx = data(3:5:end);
    vy = data(4:5:end);
    vz = data(5:5:end);
    
    vx_n = vx./sqrt(vx.^2+vy.^2);
    vy_n = vy./sqrt(vx.^2+vy.^2);
    
    vx_n(~isfinite(vx_n)) = zeros(1, sum(~isfinite(vx_n)));
    vy_n(~isfinite(vy_n)) = zeros(1, sum(~isfinite(vy_n)));
    
    
    figure
    quiver(x,y, vx_n, vy_n)
    axis([min(x), max(x), min(y), max(y)])
end