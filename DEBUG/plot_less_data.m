function [vx, vy, vz, p] = plot_less_data(trial)
       figure('units','normalized','outerposition',[0 0 1 1])
       
       subplot(1,2,1)
       plot_flow_normalized("./velocity_data/" + trial + ".txt");
       
       subplot(1,2,2)
       p = plot_pressure("./pressure_data/" + trial + ".txt");
       

end