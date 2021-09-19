function [vx, vy, vz, p] = plot_data(trial)
       figure('units','normalized','outerposition',[0 0 1 1])
       subplot(2,2,1)
       [vx,vy,vz] = plot_flow("./velocity_data/" + trial + ".txt");
       
       subplot(2,2,2)
       p = plot_pressure("./pressure_data/" + trial + ".txt");
       
       subplot(2,2,3)
       plot_flow_normalized("./velocity_data/" + trial + ".txt");
       
       subplot(2,2, 4)
       Plot_mesh("./rigid_body_data/" + trial + ".txt");
end