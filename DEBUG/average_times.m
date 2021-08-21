function [t, s, p, b, vx, vy, vz] = average_times(file_loc)
    fileID = fopen(file_loc, 'r');
    data = fscanf(fileID, '%f %f %f %f %f %f %f %f');
    fclose(fileID);
    all_s =  data(2:8:end);
    all_p =  data(3:8:end);
    all_b =  data(4:8:end);
    all_vx = data(5:8:end);
    all_vy = data(6:8:end);
    all_vz = data(7:8:end);
    all_t =  data(8:8:end);
    
    s = mean(all_s);
    p = mean(all_p);
    b = mean(all_b);
    vx = mean(all_vx);
    vy = mean(all_vy);
    vz = mean(all_vz);
    t = mean(all_t);
end