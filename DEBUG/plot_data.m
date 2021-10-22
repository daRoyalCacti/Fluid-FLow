function [vx, vy, vz, p] = plot_data(trial)
    %trial is an int
    if (trial < 10)
        trial_s = "000"+num2str(trial);
    elseif (trial < 100)
        trial_s = "00"+num2str(trial);
    elseif (trial < 1000)
        trial_s = "0"+num2str(trial);
    else
        trial_s = num2str(trial);
    end
        
    
%    figure('units','normalized','outerposition',[0 0 1 1])
   subplot(2,2,1)
   if (trial == 0)
    sgtitle({['frame = ', num2str(trial)], 't=0' })
   else 
     times_file_loc = "./times.txt";
     times_fileID = fopen(times_file_loc, 'r');
     times_data = fscanf(times_fileID, '%f %f %f %f %f %f');
     fclose(times_fileID);
     ts = times_data(1:7:end);
     sgtitle({['frame = ', num2str(trial)], ['t=', num2str(ts(trial))] })
     sgtitle({['frame = ', num2str(trial)], ['t=', num2str(ts(trial)) ]})
   end
   [vx,vy,vz] = plot_flow("./velocity_data/" + trial_s + ".txt");

   subplot(2,2,2)
   p = plot_pressure("./pressure_data/" + trial_s + ".txt");

   subplot(2,2,3)
%    plot_flow_normalized("./velocity_data/" + trial_s + ".txt");
   plot_flow_stream("./velocity_data/" + trial_s + ".txt", "./pressure_data/" + trial_s + ".txt");

   subplot(2,2, 4)
   Plot_mesh("./rigid_body_data/" + trial_s + ".txt");
end