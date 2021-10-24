no_files = length(dirPlus('velocity_data'));

figure('units','normalized','outerposition',[0 0 1 1]);
axis tight manual
%set(gcf, 'Position', get(0, 'Screensize'));
hold on
v = VideoWriter('test.avi', 'Motion JPEG AVI');
v.Quality = 40; %lower quality to save on disk space
v.FrameRate = 60;
open(v)

for frame = 0:no_files-1
    clf;
    plot_data(frame);
    framed = getframe(gcf);
    writeVideo(v,framed);
end

close(v)
close(gcf)