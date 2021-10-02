function [vx, vy, vz, p] = plot_less_data(trial)
%        figure('units','normalized','outerposition',[0 0 1 1])
        figure('Position',[0 0 1000 500])
       
       subplot(1,2,1)
       plot_flow_norm_less("./velocity_data/" + trial + ".txt");
       
       subplot(1,2,2)
       p = plot_pressure("./pressure_data/" + trial + ".txt");
       

end