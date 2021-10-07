fprintf("Finding files\n")
All_v_files = dirPlus('velocity_data');   %https://stackoverflow.com/questions/2652630/how-to-get-all-files-under-a-specific-directory-in-matlab#2654459
All_p_files = dirPlus('pressure_data'); 
All_b_files = dirPlus('rigid_body_data');

no_files = length(All_v_files);

no_files = 174; %for testing

fprintf("Reading files\n")
data = struct;
for ii = 1:no_files
    [data(ii).x, data(ii).y , data(ii).vx , data(ii).vy , data(ii).vz ] = read_v_data( cell2mat(All_v_files(ii)) );
    
    [data(ii).p, data(ii).dx, data(ii).dy, data(ii).x_sort, data(ii).y_sort] = read_p_data( cell2mat(All_p_files(ii)) );
    [data(ii).bx, data(ii).by, data(ii).bz] = read_b_data( cell2mat(All_b_files(ii)) );
end


fprintf("Finding mins and maxs\n")
minx = min(data(1).x);
miny = min(data(1).y);

maxx = max(data(1).x);
maxy = max(data(1).y);

max_vx = -inf;
max_vy = -inf;

min_p = inf;
max_p = -inf;

for ii = 1:no_files
    max_vx = max([data(ii).vx; max_vx]);
    max_vy = max([data(ii).vy; max_vy]);
    
    min_p = min([min(data(ii).p), min_p]);
    max_p = max([max(data(ii).p), max_p]);
end


figure('units','normalized','outerposition',[0 0 1 1]);
axis tight manual
%set(gcf, 'Position', get(0, 'Screensize'));
hold on

%creating video writer
%creating a .avi file with 'motion jpeg'
% - preffered would be an mp4 with h264 but this is not
% available on linux
v = VideoWriter('test.avi', 'Motion JPEG AVI');
v.Quality = 40; %lower quality to save on disk space
v.FrameRate = 30;
open(v)

for frame = 1:no_files
    fprintf("Frame %i/%i\n", frame, no_files);
    %fprintf('drawing frame %d of %d\n', frame+1, frames)
    
    clf
    sgtitle(['frame = ', num2str(frame)])
    subplot(2,2,1)
    quiver(data(frame).x, data(frame).y, data(frame).vx/max_vx, data(frame).vy/max_vy)
    axis([minx, maxx, miny, maxy])
    xlabel('x')
    ylabel('y')

    
    subplot(2,2,2)
    h = heatmap(unique(x),unique(y), flip(data(frame).p), 'GridVisible', 'off');
    
    %changing the axes
    Ax = gca;
    
    xlabels = nan(size(Ax.XDisplayData));
    xlabels(1:10:end) = floor( linspace(data(frame).x_sort(1), data(frame).x_sort(end), length(test(1:10:end)) ) *100 )/100;
    Ax.XDisplayLabels = xlabels;
    
    ylabels = nan(size(Ax.YDisplayData));
    ylabels(end:-10:1) = floor( linspace(data(frame).y_sort(1), data(frame).y_sort(end), length(test(1:10:end)) ) *100 )/100;
    Ax.YDisplayLabels = ylabels;
    
    Ax.XLabel = 'x';
    Ax.YLabel = 'y';
    Ax.MissingDataLabel = '';
    
    %rotating x tick labels
    set(struct(h).NodeChildren(3), 'XTickLabelRotation', 0);
     
     
    
    subplot(2,2,3)
    quiver(data(frame).x, data(frame).y, data(frame).vx./sqrt( data(frame).vx.^2 + data(frame).vy.^2), data(frame).vy./sqrt( data(frame).vx.^2 + data(frame).vy.^2))
%     starty = linspace(min(data(frame).y), max(data(frame).y), 20);
%     startx = ones(size(starty)) * min(data(frame).x);
%     streamline(data(frame).x, data(frame).y, data(frame).vx, data(frame).vy, startx, starty)
    axis([minx, maxx, miny, maxy])
    xlabel('x')
    ylabel('y')
    
    subplot(2,2,4)
    plot3(data(frame).bx, data(frame).by, data(frame).bz, 'o')
    xlabel('x')
    ylabel('y')
    zlabel('z')


    % Capture the plot as an image
    framed = getframe(gcf);
    writeVideo(v,framed);
end

close(v)
close(gcf)


function [x,y, vx, vy, vz] = read_v_data(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f');
    fclose(fileID);
    x = data(1:5:end);
    y = data(2:5:end);
    
    vx = data(3:5:end);
    vy = data(4:5:end);
    vz = data(5:5:end);

end


function [p, dx, dy, x_sort, y_sort] = read_p_data(file_loc)
fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f');
    fclose(fileID);
    x = data(1:3:end);
    y = data(2:3:end); 
    p = data(3:3:end);
        
    x_sort = sort(unique(x));
    dx = x_sort(2)-x_sort(1);
    y_sort = sort(unique(y));
    dy = y_sort(2)-y_sort(1);
    
        off = 0;
    p_r = p;
    for ii = 1:length(p)-1
        ind = ii+off;
        if (x(ii+1) > x(ii))
            if (x(ii+1) - x(ii) > 1.5*dx)
                no_points = floor( (x(ii+1) - x(ii))/dx );
                p_r = [p_r(1:ind); nan*(1:no_points)'; p_r(ind+1:end) ];
                off = off + no_points;
            end
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

end


function [x,y,z] = read_b_data(file_loc)
fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f');
    fclose(fileID);
    x = data(1:3:end);
    y = data(2:3:end);
    z = data(3:3:end);
end


