function [vx, vy, vz, p] = plot_data(trial)
       figure
       subplot(2,1,1)
       [vx,vy,vz] = plot_flow("./velocity_data/" + trial + ".txt");
       subplot(2,1,2)
       
       p = plot_pressure("./pressure_data/" + trial + ".txt");
end