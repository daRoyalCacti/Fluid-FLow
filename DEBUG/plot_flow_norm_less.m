%plots the flow for a give time
function [vx, vy, vz,x,y] = plot_flow_norm_less(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f');
    fclose(fileID);
    x = data(1:5:end);
    y = data(2:5:end);
    
    vx = data(3:5:end);
    vy = data(4:5:end);
    vz = data(5:5:end);
    
    skip = 6;
    
    ux = unique(x);
    uy = unique(y);
    
    xs = ux(1:skip:end);
    ys = uy(1:skip:end);
    
    inds_x = zeros(size(x));
    inds_y = zeros(size(x));
    for ii = 1:length(x)
        inds_x(ii) = ismember(x(ii),xs);
        inds_y(ii) = ismember(y(ii),ys);
    end
 
    
    inds = inds_x & inds_y;
    
    x = x(inds);
    y = y(inds);
    vx = vx(inds);
    vy = vy(inds);
    
    vx_n = vx./sqrt(vx.^2+vy.^2);
    vy_n = vy./sqrt(vx.^2+vy.^2);
    
    vx_n(~isfinite(vx_n)) = zeros(1, sum(~isfinite(vx_n)));
    vy_n(~isfinite(vy_n)) = zeros(1, sum(~isfinite(vy_n)));
    
    
%     figure
    quiver(x,y, vx_n, vy_n)
    axis([min(x), max(x), min(y), max(y)])
    xlabel('x')
    ylabel('y')
end