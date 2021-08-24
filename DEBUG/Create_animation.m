All_files = dirPlus('rigid_body_data');   %https://stackoverflow.com/questions/2652630/how-to-get-all-files-under-a-specific-directory-in-matlab#2654459

data = struct;
for ii = 1:length(All_files)
    [x,y,z] = read_data( cell2mat(All_files(ii)) );
    data(ii).x = x;
    data(ii).y = y;
    data(ii).z = z;
end

minx = inf;
miny = inf;
minz = inf;

maxx = -inf;
maxy = -inf;
maxz = -inf;

for ii = 1:length(All_files)
    minx = min([data(ii).x; minx]);
    miny = min([data(ii).y; miny]);
    minz = min([data(ii).z; minz]);
    
    maxx = max([data(ii).x; maxx]);
    maxy = max([data(ii).y; maxy]);
    maxz = max([data(ii).z; maxz]);
end

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
open(v)

for frame = 1:length(All_files)
    %fprintf('drawing frame %d of %d\n', frame+1, frames)
    
    clf
    %plot3(x(:,frame+1),y(:,frame+1),z(:,frame+1), 'ro', 'MarkerSize',3)
    plot3(data(frame).x, data(frame).y, data(frame).z, 'o')
    view(-45, 45)
    axis([minx, maxx, miny, maxy, minz, maxz])
    xlabel('x')
    ylabel('y')
    zlabel('z')


    % Capture the plot as an image
    framed = getframe(gcf);
    writeVideo(v,framed);
end

close(v)
close(gcf)


function [x,y,z] = read_data(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f');
    fclose(fileID);
    x = data(1:3:end);
    y = data(2:3:end);
    z = data(3:3:end);
end
