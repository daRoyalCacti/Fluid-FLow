All_files = dirPlus('plotting_raw_data');   %https://stackoverflow.com/questions/2652630/how-to-get-all-files-under-a-specific-directory-in-matlab#2654459

data = struct;
for ii = 1:length(All_files)
    [x,y,vx, vy, vz] = read_data( cell2mat(All_files(ii)) );
    data(ii).x = x;
    data(ii).y = y;
    data(ii).vx = vx;
    data(ii).vy = vy;
    data(ii).vz = vz;
end

minx = min(x);
miny = min(y);

maxx = max(x);
maxy = max(y);

max_vx = -inf;
max_vy = -inf;

for ii = 1:length(All_files)
    max_vx = max([data(ii).vx; max_vx]);
    max_vy = max([data(ii).vy; max_vy]);
end
% max_vx = 1;
% max_vy = 1;

figure;
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

for frame = 1:length(All_files)
    %fprintf('drawing frame %d of %d\n', frame+1, frames)
    
    clf
    %plot3(x(:,frame+1),y(:,frame+1),z(:,frame+1), 'ro', 'MarkerSize',3)
    %plot3(data(frame).x, data(frame).y, data(frame).z, 'o')
    quiver(data(frame).x, data(frame).y, data(frame).vx/max_vx, data(frame).vy/max_vy)
    axis([minx, maxx, miny, maxy])
    xlabel('x')
    ylabel('y')


    % Capture the plot as an image
    framed = getframe(gcf);
    writeVideo(v,framed);
end

close(v)
close(gcf)

%plots the flow for a give time
function [x,y, vx, vy, vz] = read_data(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f');
    fclose(fileID);
    x = data(1:5:end);
    y = data(2:5:end);
    
    vx = data(3:5:end);
    vy = data(4:5:end);
    vz = data(5:5:end);

end


