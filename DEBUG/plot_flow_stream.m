%plots the flow for a give time
function [vx, vy, vz] = plot_flow_stream(file_loc, pressure_file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f');
    fclose(fileID);
%     x = data(1:5:end);
%     y = data(2:5:end);
    
    vx = data(3:5:end);
    vy = data(4:5:end);
    vz = data(5:5:end);
    
    
    fileIDp = fopen(pressure_file_loc, 'r');
    datap = fscanf(fileIDp, '%f %f %f');
    fclose(fileIDp);
    x = datap(1:3:end);
    y = datap(2:3:end); 
    
    
    
    
    x_sort = sort(unique(x));
    dx = x_sort(2)-x_sort(1);
    y_sort = sort(unique(y));
    dy = y_sort(2)-y_sort(1);
    
    off = 0;
    
    x = [x_sort(1); x; x_sort(end)];
    y=[y_sort(1); y; y_sort(end)];
    x_r = x;
    y_r = y;
    
    
    
    vx = [0; vx; 0];
    vy = [0; vy; 0];
    len = length(vx)-1;
    for ii = 1:len
        ind = ii+off;
        if (x(ii+1) - x(ii) > 1.5*dx)
            no_points = round( (x(ii+1) - x(ii))/dx )-1;%length(x(ii):dx:x(ii+1))-2;%floor( (x(ii+1) - x(ii))/dx );
            
            x_repl = linspace(x(ii)+dx, x(ii+1)-dx, no_points);%(x(ii)):dx:(x(ii+1)-dx);
            x_repl = x_repl';
            x_r = [x_r(1:ind); x_repl; x_r(ind+1:end) ];
            y_r = [y_r(1:ind); y_r(ind)*ones(size(x_repl)); y_r(ind+1:end)];
            
            vx = [vx(1:ind); 0*(1:no_points)'; vx(ind+1:end) ];
            vy = [vy(1:ind); 0*(1:no_points)'; vy(ind+1:end) ];
            off = off + no_points;
        end
        
        ind = ii+off;
        if ( (y(ii+1) > y(ii)) && abs( x(ii+1) - x_sort(1) ) > dx/3  )
            x_r = [x_r(1:ind); x_sort(1); x_r(ind+1:end)];
            y_r = [y_r(1:ind); y(ii+1); y_r(ind+1:end)];
            
            vx = [vx(1:ind); 0; vx(ind+1:end) ];
            vy = [vy(1:ind); 0; vy(ind+1:end) ];
            
            off = off+1;
        end
        
        ind = ii+off;
        if ( (y(ii+1) > y(ii)) && abs( x(ii) - x_sort(end) ) > dx/3  )
            x_r = [x_r(1:ind); x_sort(end); x_r(ind+1:end)];
            y_r = [y_r(1:ind); y(ii); y_r(ind+1:end)];
            
            vx = [vx(1:ind); 0; vx(ind+1:end) ];
            vy = [vy(1:ind); 0; vy(ind+1:end) ];
            
            off = off+1;
        end
    end
    
    w = sqrt(length(vx));


    xp = reshape(x_r, w, w)';
    yp = reshape(y_r, w, w)';
    vxp = reshape(vx, w, w)';
    vyp = reshape(vy, w, w)';
    
    
    streamslice(xp, yp, vxp, vyp, 3, 'cubic', 'arrows')
    axis([x_sort(1), x_sort(end), y_sort(1), y_sort(end)])
















