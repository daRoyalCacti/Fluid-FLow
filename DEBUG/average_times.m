% function average_times      %function so workspace doesn't get filled with variables

file_loc = "./times.txt";
fileID = fopen(file_loc, 'r');
data = fscanf(fileID, '%f %f %f %f %f %f');
fclose(fileID);
all_s =  data(2:7:end);
all_p =  data(3:7:end);
all_b =  data(4:7:end);
all_v = data(5:7:end);
all_mesh =  data(6:7:end);
all_t =  data(7:7:end);

s = mean(all_s);
p = mean(all_p);
b = mean(all_b);
v = mean(all_v);
t = mean(all_t);
m = mean(all_mesh);

pie([m, s, p, b, v],[1,0,0,0,0], {'updating mesh', 'making s', 'solving p', 'making b', 'solving v'})
text(-0.5,-1.3, ['average time = ', num2str(t), 's'])
