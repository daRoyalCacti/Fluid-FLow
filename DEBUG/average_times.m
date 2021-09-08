function average_times      %function so workspace doesn't get filled with variables

file_loc = "./times.txt";
fileID = fopen(file_loc, 'r');
data = fscanf(fileID, '%f %f %f %f %f %f %f %f %f');
fclose(fileID);
all_s =  data(2:9:end);
all_p =  data(3:9:end);
all_b =  data(4:9:end);
all_vx = data(5:9:end);
all_vy = data(6:9:end);
all_vz = data(7:9:end);
all_mesh =  data(8:9:end);
all_t =  data(9:9:end);

s = mean(all_s);
p = mean(all_p);
b = mean(all_b);
vx = mean(all_vx);
vy = mean(all_vy);
vz = mean(all_vz);
t = mean(all_t);
m = mean(all_mesh);

pie([m, s, p, b, vx, vy, vz],[1,0,0,0,0,0,0], {'updating mesh', 'making s', 'solving p', 'making b', 'solving vx', 'solving vy', 'solving vz'})
text(-0.5,-1.3, ['average time = ', num2str(t), 's'])
