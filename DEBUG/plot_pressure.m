function p = plot_pressure(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f');
    fclose(fileID);
    x = data(1:3:end);
    y = data(2:3:end); 
    p = data(3:3:end);
    
%     axis([min(x), max(x), min(y), max(y)])
    
    x_sort = sort(unique(x));
    dx = x_sort(2)-x_sort(1);
    y_sort = sort(unique(y));
    dy = y_sort(2)-y_sort(1);
   
%     if ( max(p) - min(p) > 0.00000001)
%         p_n = ( p - min(p) ) / (max(p) - min(p) );
%     else
%         p_n = zeros(size(p));
%     end
%             
%     interp_col1 = [100, 100, 100]/100;
%     interp_col2 = [3.5, 43.9, 77.3]/100;
%     
%     interp_func = @(x) interp_col1*x + interp_col2*(1-x);
% 
%     
%     for ii = 1:(length(x))
% %         rectangle('Position',[x(ii),y(ii),dx,dy],'FaceColor',[0 0 p_n(ii)],'EdgeColor',[0 0 p_n(ii)],'LineWidth',0.001)
%         rectangle('Position',[x(ii),y(ii),dx,dy],'FaceColor',interp_func(p_n(ii)),'EdgeColor', interp_func(p_n(ii)), 'LineWidth',0.001)
%     end
%     
%     xlabel('x')
%     ylabel('y')
%     c=colorbar;
%     c.Limits = [0, 1];
%     %setting the colour bar colours
%     map = [linspace(interp_col1(1) ,interp_col2(1) , 100)', linspace(interp_col1(2) ,interp_col2(2) , 100)', linspace(interp_col1(3) ,interp_col2(3) , 100)'];
%     colormap(map)
%     
%     %labelling the colour bar
%     no_ticks = length(c.TickLabels);
%     label_vals = linspace(min(p), max(p), no_ticks);
%     label_vals = round(label_vals*100)/100; %giving the labels only 2 decimal places
%     for ii = 1:no_ticks
%         tick_labels(ii) = {num2str(label_vals(ii))};
%     end
%     c.TickLabels = tick_labels;
%     
%      rectangle('Position',[min(x),min(y),max(x)-min(x),max(y)-min(y)],'EdgeColor', [0,0,0], 'LineWidth',1.5)
    
    off = 1;
    p_r = [nan; p; nan];
    for ii = 1:length(p)-1
        ind = ii+off;
        if (x(ii+1) - x(ii) > 1.5*dx)
            no_points = length(x(ii):dx:x(ii+1))-2;%floor( (x(ii+1) - x(ii))/dx );
            p_r = [p_r(1:ind); nan*(1:no_points)'; p_r(ind+1:end) ];
            off = off + no_points;
        end
        
        ind = ii+off;
        if ( (y(ii+1) > y(ii)) && abs( x(ii+1) - x_sort(1) ) > dx/3  )
            p_r = [p_r(1:ind); nan; p_r(ind+1:end)];
            off = off+1;
        end
        
        ind = ii+off;
        if ( (y(ii+1) > y(ii)) && abs( x(ii) - x_sort(end) ) > dx/3  )
            p_r = [p_r(1:ind); nan; p_r(ind+1:end)];
            off = off+1;
        end
        
    end
    
    l = length(p_r);
    if ( sqrt(l) ~= floor(sqrt(l)))
        error("assumption that grid is a square failed")
    end
    
    p = zeros(sqrt(l), sqrt(l));
    
    for ii = 1:sqrt(l)
        for jj = 1:sqrt(l)
            p(jj,ii) = p_r(ii + sqrt(l)*(jj-1) );
        end
    end
     
    
    
%     figure
    h = heatmap(unique(x),unique(y), flip(p), 'GridVisible', 'off');
%     h = heatmap(unique(y),unique(x), flip(p'), 'GridVisible', 'off');
    
    %changing the axes
    Ax = gca;
    
    xlabels = nan(size(Ax.XDisplayData));
    xlabels(1:10:end) = floor( linspace(x_sort(1), x_sort(end), length(xlabels(1:10:end)) ) *100 )/100;
    Ax.XDisplayLabels = xlabels;
    
    ylabels = nan(size(Ax.YDisplayData));
    ylabels(end:-10:1) = floor( linspace(y_sort(1), y_sort(end), length(ylabels(1:10:end)) ) *100 )/100;
    Ax.YDisplayLabels = ylabels;
    
    Ax.XLabel = 'x';
    Ax.YLabel = 'y';
    Ax.MissingDataLabel = '';
    
    %rotating x tick labels
    set(struct(h).NodeChildren(3), 'XTickLabelRotation', 0);
    
%     figure
%     [x_,y_] = meshgrid(unique(x), unique(y) );
%     mesh(x_, y_, p)
%     xlabel('x')
%     ylabel('y')
%     zlabel('p')
% end