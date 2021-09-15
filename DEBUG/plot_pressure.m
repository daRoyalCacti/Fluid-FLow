function p_r = plot_pressure(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f');
    fclose(fileID);
    x = data(1:3:end);
    y = data(2:3:end); 
    p_r = data(3:3:end);
    
    axis([min(x), max(x), min(y), max(y)])
    
    x_sort = sort(unique(x));
    dx = x_sort(2)-x_sort(1);
    y_sort = sort(unique(y));
    dy = y_sort(2)-y_sort(1);
    
    p_n = ( p_r - min(p_r) ) / (max(p_r) - min(p_r) );
    
    
    for ii = 1:(length(x)-1)
        rectangle('Position',[x(ii),y(ii),dx,dy],'FaceColor',[0 0 p_n(ii)],'EdgeColor',[0 0 p_n(ii)],'LineWidth',0.001)
    end

%     l = length(p_r);
%     if ( sqrt(l) ~= floor(sqrt(l)))
%         error("assumption that grid is a square failed")
%     end
%     
%     p = zeros(sqrt(l), sqrt(l));
%     
%     for ii = 1:sqrt(l)
%         for jj = 1:sqrt(l)
%             p(ii,jj) = p_r(ii + sqrt(l)*(jj-1) );
%         end
%     end
%      
%     
%     
% %     figure
%     heatmap(unique(x),unique(y), flip(p), 'GridVisible', 'off');
%     
%     %turning off axes
%     Ax = gca;
%     Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
%     Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
end