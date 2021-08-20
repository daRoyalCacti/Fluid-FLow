fprintf("Finding files\n")
All_v_files = dirPlus('velocity_data');   %https://stackoverflow.com/questions/2652630/how-to-get-all-files-under-a-specific-directory-in-matlab#2654459
All_p_files = dirPlus('pressure_data'); 

no_files = length(All_v_files);

fprintf("Reading files\n")
data = struct;
for ii = 1:no_files
    [data(ii).x, data(ii).y , data(ii).vx , data(ii).vy , data(ii).vz ] = read_v_data( cell2mat(All_v_files(ii)) );
    
    data(ii).p = read_p_data( cell2mat(All_p_files(ii)) );
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
    subplot(2,1,1)
    quiver(data(frame).x, data(frame).y, data(frame).vx/max_vx, data(frame).vy/max_vy)
    axis([minx, maxx, miny, maxy])
    xlabel('x')
    ylabel('y')
    
    subplot(2,1,2)
    heatmap(unique(data(frame).x),unique(data(frame).y), flip(data(frame).p), 'GridVisible', 'off', 'ColorLimits', [min_p, max_p]);
    
    %turning off axes
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));


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


function p = read_p_data(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f');
    fclose(fileID);
    %x = data(1:3:end);
    %y = data(2:3:end); 
    p_r = data(3:3:end);
    
%     p = data;
    l = length(p_r);
    if ( sqrt(l) ~= floor(sqrt(l)))
        error("assumption that grid is a square failed")
    end
    
    p = zeros(sqrt(l), sqrt(l));
    
    for ii = 1:sqrt(l)
        for jj = 1:sqrt(l)
            p(ii,jj) = p_r(ii + sqrt(l)*(jj-1) );
        end
    end

end


