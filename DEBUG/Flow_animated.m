fprintf("Finding files\n")
All_v_files = dirPlus('velocity_data');   %https://stackoverflow.com/questions/2652630/how-to-get-all-files-under-a-specific-directory-in-matlab#2654459
All_p_files = dirPlus('pressure_data'); 
All_b_files = dirPlus('rigid_body_data');

% no_files = length(All_v_files);
fprintf('no_files has been changed for testing\n')
no_files = 218;

fprintf("Reading files\n")
data = struct;
for ii = 1:no_files
    [data(ii).x, data(ii).y , data(ii).vx , data(ii).vy , data(ii).vz ] = read_v_data( cell2mat(All_v_files(ii)) );
    
    [data(ii).p, data(ii).dx, data(ii).dy] = read_p_data( cell2mat(All_p_files(ii)) );
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
    for ii = 1:(length(x))
        rectangle('Position',[data(frame).x(ii),data(frame).y(ii),data(frame).dx,data(frame).dy],'FaceColor',...
            [0 data(frame).p(ii) 0],'EdgeColor',[0 data(frame).p(ii) 0],'LineWidth',0.001)
    end
    
    subplot(2,2,3)
    quiver(data(frame).x, data(frame).y, data(frame).vx./sqrt( data(frame).vx.^2 + data(frame).vy.^2), data(frame).vy./sqrt( data(frame).vx.^2 + data(frame).vy.^2))
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


function [p_n, dx, dy] = read_p_data(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f');
    fclose(fileID);
    x = data(1:3:end);
    y = data(2:3:end); 
    p = data(3:3:end);
    
    axis([min(x), max(x), min(y), max(y)])
    
    x_sort = sort(unique(x));
    dx = x_sort(2)-x_sort(1);
    y_sort = sort(unique(y));
    dy = y_sort(2)-y_sort(1);
   
    if ( max(p) - min(p) > 0.00000001)
        p_n = ( p - min(p) ) / (max(p) - min(p) );
    else
        p_n = zeros(size(p));
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


