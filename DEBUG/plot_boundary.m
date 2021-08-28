%plots the flow for a give time
function [vx, vy, vz] = plot_boundary(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f');
    fclose(fileID);
    x = data(1:3:end);
    y = data(2:3:end);
    
    b = data(3:3:end);
    
    l = length(b);
    if ( sqrt(l) ~= floor(sqrt(l)))
        error("assumption that grid is a square failed")
    end
    
    p = zeros(sqrt(l), sqrt(l));
    
    for ii = 1:sqrt(l)
        for jj = 1:sqrt(l)
            p(ii,jj) = b(ii + sqrt(l)*(jj-1) );
        end
    end
     
    
    
%     figure
    heatmap(unique(x),unique(y), flip(p), 'GridVisible', 'off');
    
    %turning off axes
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
end