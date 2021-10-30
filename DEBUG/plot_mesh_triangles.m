function [i,v,c] = plot_mesh_triangles(file_loc, output_loc)

vert_file = "vert_data.txt";
ind_file = "ind_data.txt";

fileID = fopen(vert_file, 'r');
vr = fscanf(fileID, '%f %f %f');
fclose(fileID);
v = [vr(1:3:end), vr(2:3:end), vr(3:3:end) ];

fileID = fopen(ind_file, 'r');
ir = fscanf(fileID, '%f %f %f');
fclose(fileID);
i = [ir(1:3:end), ir(2:3:end), ir(3:3:end) ] +1;




%Read in the pressure data
% file_loc = "pressure_dump/0001.txt";
fileID = fopen(file_loc, 'r');
pressure_data_raw = fscanf(fileID, '%f %f %f %f');
fclose(fileID);

px = pressure_data_raw(1:4:end);
py = pressure_data_raw(2:4:end);
pz = pressure_data_raw(3:4:end);
p = pressure_data_raw(4:4:end);

clear pressure_data_raw

ls = max(v) - min(v);
ml = max(ls);

ls2 = ls;
ls2( ls2 < ml/2 ) = ml/2;

ds = 2*ls2 / 129; 
dx = ds(1);
dy = ds(2);
dz = ds(3);



c = zeros(length(v), 1);

%find values to draw pressure at
for ii = 1:length(v)
    fprintf('%i/%i\n', ii, length(v))
        vec = v(ii, :);
        %find pressure at inds
        inds1 = (px > vec(1)) & (px < vec(1)+dx);
        inds2 = (py > vec(2)) & (py < vec(2)+dy);
        inds3 = (pz > vec(3)) & (pz < vec(3)+dz);

        indst = inds1 & inds2 & inds3;
        
        if (sum(indst) == 0)
            for inc = 2:1000
                if (sum(inds1) == 0)
                    inds1 = (px > vec(1)) & (px < vec(1)+2*dx);
                end
                if (sum(inds2) == 0)
                    inds2 = (py > vec(2)) & (py < vec(2)+2*dy);
                end
                if (sums(inds3) == 0)
                    inds3 = (pz > vec(3)) & (pz < vec(3)+2*dz);
                end

                indst = inds1 & inds2 & inds3;
                
                if (sum(indst) > 0)
                    break;
                end
            
            end
            
            c(ii) = mean(p(indst));
        else 
            c(ii) = mean(p(indst)); %mean just in case indst has 2 values in it
        end
        
        
end






figure;
% drawMesh(v, i, c)
colormap hot
patch('Faces', i, 'Vertices', v, 'FaceVertexCData',c, 'FaceColor','interp', 'EdgeAlpha', 0.5)
view(3)
daspect([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')
colorbar
view([180+45 45])

savefig(gcf, [output_loc, '/raw_press_dist.fig'])
exportgraphics(gcf, [output_loc, '/press_dist.png'])


end






