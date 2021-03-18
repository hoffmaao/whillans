import numpy as np
import pandas as pd
from math import pi
from constants import (cice, cair)

def RollingRadon(data_x_or_filename,data_y,Data,window,angle_thresh,
    plotter,surface_bottom,movie_flag,max_frequency,o_f_vertical=6,o_f_horizontal = 2,snr_thresh = 2,snr_fac = 1,vr = 3):

if isstr(data_x_or_filename) == 1
    
    load(data_x_or_filename)
    
    if exist('x') == 0
        [x y] = polarstereo_fwd(Latitude,Longitude);
    end
    if Time(2)-Time(1) > 1e-6
        Time = Time*10^-6;
    end
    
    dist = distance_vector(x,y);
    data_y = Time;
    
    Data = lp(Data);
    
else
    dist = data_x_or_filename; 
    
    if surface_bottom ~= 0
        surf_bot_dim = size(surface_bottom);
        if min(surf_bot_dim) == 1
            Surface = surface_bottom;
        else
            if surf_bot_dim(1) == min(surf_bot_dim)
                Surface = surface_bottom(1,:);
                Bottom = surface_bottom(2,:);
            else
                Surface = surface_bottom(:,1);
                Bottom = surface_bottom(:,2);
            end
        end
    end
end

%% Set-up the loop
filt_data = lp(Data);

clearvars -except Data dist data_y Surface Bottom window_size ...
    window_size2 o_f_vertical o_f_horizontal snr_thresh plotter ...
    movie_flag snr_fac max_frequency xstep_roll ystep_roll ...
    angle_thresh pr vr
;



overload_factor = 1000; 

if length(Data(1,:)) > overload_factor
  steps = ceil(length(Data(1,:))/overload_factor);
  breaks = [1:overload_factor:(length(Data(1,:))+1) (length(Data(1,:))+1)];
else
    steps = 1;
    breaks = [1 length(Data(1,:))];
end


    slope_colors = b2r2(-angle_thresh(1),angle_thresh(2));
    slope_vals = -angle_thresh(1): ...
        (2*angle_thresh(1)+1)/length(slope_colors(:,1)):angle_thresh(1);
    total_time = 0;
    
for k = 1:steps
    %%%% Update to the console
    disp(['Starting Window ',num2str(k),' of ',num2str(steps), ...
        ', Total Time - ',sprintf('%.02f',total_time),' min'])
    
    %% Data Preconditioning:
    clearvars xaxis yaxis data
    

    filt_data = Data(:,breaks(k):breaks(k+1)-1);
    if exist('max_frequency') == 1
        [data xaxis yaxis] = regrid(dist(breaks(k):breaks(k+1)-1), ...
            data_y,filt_data,1,max_frequency);
        time = yaxis;
        yaxis = yaxis*cice/2;
    else
        [data xaxis yaxis] = regrid(dist(breaks(k):breaks(k+1)-1), ...
            data_y,filt_data,0,0);
    end
    
    if exist('Bottom') == 1
        Bottom2 = interp1(dist(breaks(k):breaks(k+1)-1), ...
            interpNaN(Bottom(breaks(k):breaks(k+1)-1)),xaxis);
        Surface2 = interp1(dist(breaks(k):breaks(k+1)-1), ...
            interpNaN(Surface(breaks(k):breaks(k+1)-1)),xaxis);
    end
    
    if exist('previous_xsteps') == 0
        previous_xsteps = 0;
    end
    
    if k < steps
        roll_steps = round((length(data(1,:))-window_size)/xstep_roll);
    else
        roll_steps = floor((length(data(1,:))-window_size)/xstep_roll);
    end

    roll_steps2 = round((length(data(:,1))-window_size2)/ystep_roll);

    keep_val = 1;
    
    if exist('opt_angle') == 0
        opt_angle = zeros(roll_steps2,roll_steps)*NaN;
        status_flag = zeros(size(opt_angle));
        means = zeros(size(opt_angle));
    else
        opt_angle = [opt_angle zeros(roll_steps2,roll_steps)*NaN];
        status_flag = [status_flag zeros(roll_steps2,roll_steps)];
        means = [means zeros(roll_steps2,roll_steps)];
    end
    
    
    end_counter1 = 0;
    end_counter2 = 0;
    snr_counter = 1;
    counter1 = 1;
    updater = 10;
    
    tic
    
    for i = 1:roll_steps
        
        %%% Restart the variability record
        variability_record = [];
        last_val = [];
        
        %%% Determine the Window area for the horizontal dimension
        start = (i-1)*xstep_roll+1;
        stop = min([start+window_size-1 length(data(1,:))]);
        opt_x(i+previous_xsteps) = xaxis(stop-floor((window_size-1)/2));
        if i > 1
            if opt_x(end) <= opt_x(end-1)
            end_counter1 = end_counter1+1;
            opt_x(end) = xaxis(stop-floor((window_size-1)/2)+end_counter1);
            else
            end_counter1 = 0;
            end
        end
        
        centerx_ind = stop-floor((window_size-1)/2);
        
        power_dist = conv(data(:,centerx_ind), ...
            ones(round(length(data(:,centerx_ind))/50),1),'same')./ ...
            conv(ones(size(data(:,centerx_ind))), ...
            ones(round(length(data(:,centerx_ind))/50),1),'same');
        power_dist_mean = mean(power_dist);
        power_dist_std = std(power_dist);

        
        for j = 1:roll_steps2
            start2 = (j-1)*floor(window_size2/o_f_vertical)+1;
            stop2 = min([(j-1)*floor(window_size2/o_f_vertical) + ...
                window_size2 length(data(:,1))]);
            opt_y(j) = yaxis(stop2-floor((stop2-start2)/2));
            centery_ind = stop2-floor((stop2-start2)/2);
            
      if i > 1
      if opt_x(end) <= opt_x(end-1)
          end_counter2 = end_counter2+1;
          opt_x(end) = xaxis(stop2-floor((window_size2-1)/2)+end_counter2);
      else
          end_counter2 = 0;
      end
      end
            
            if exist('Bottom') == 1
                if time(centery_ind) < Bottom2(centerx_ind) & ...
                        time(centery_ind) > Surface2(centerx_ind)
                    skipflag = 0;
                else
                    skipflag = 1;
                    status_flag(j,i+previous_xsteps) = 1;
                end
            else
                skipflag = 0;
            end
            
        if skipflag == 0;

            radon_data = data(start2:stop2,start:stop);
            means(j,i+previous_xsteps) = mean(mean(radon_data));

            snr_win_data = data(start2:stop2,centerx_ind);
            snr_std = std(snr_win_data);
            snr = 2*snr_std;

            if snr < snr_thresh | power_dist(centery_ind) < ...
                    power_dist_mean - snr_fac*power_dist_std
                opt_angle(j,i+previous_xsteps) = NaN;
                skipflag = 1;
                status_flag(j,i+previous_xsteps) = 2;
            else

                if length(xaxis(start:stop)) == 0
                    keyboard
                end

                %% Compute the Radon transform
                [opt_angle(j,i+previous_xsteps) rd trash trash rsnr] = ...
                    radon_ndh(xaxis(start:stop),yaxis(start2:stop2), ...
                    radon_data,angle_thresh(1),0,0);
                if isnan(opt_angle(j,i+previous_xsteps)) == 1
                    status_flag(j,i+previous_xsteps) = 2;
                end

                %%% This identifies if the value exceeds the second
                %%% entry in angle_thresh
                if abs(opt_angle(j,i+previous_xsteps)) > angle_thresh(2)
                    status_flag(j,i+previous_xsteps) = 3;
                    opt_angle(j,i+previous_xsteps) = NaN;
                end                        
            end
        else
                opt_angle(j,i+previous_xsteps) = NaN;
        end
            
            
    if skipflag == 0
        variability_thresh = 4;
                if j~= 1 && length(last_val) > 0     
                    if abs(last_val-opt_angle(j,i+previous_xsteps)) > ...
                            variability_thresh
                        variability_record = [variability_recor
                            d ...
                            opt_angle(j,i+previous_xsteps)];
                        if length(variability_record) >= vr; 
                            opt_angle(j-vr+1:j,i+previous_xsteps) = ...
                                variability_record(1:vr);
                            last_val = opt_angle(j,i+previous_xsteps);
                            status_flag(j-vr+1:j,i+previous_xsteps) = 3;
                        else                                
                            opt_angle(j,i+previous_xsteps) = last_val;
                            last_val = opt_angle(j,i+previous_xsteps);
                            status_flag(j,i+previous_xsteps) = 3;
                        end
                    else
                        variability_record = [];
                    end
                else
                    last_val = opt_angle(j,i+previous_xsteps);
                    variability_record = [];
                end
    end
end




%% Produce the final results image
zero_inds = find(opt_x ~= 0);
slope_x = opt_x(zero_inds);
slope_y = opt_y;
slopes = opt_angle(:,zero_inds);
means = means(:,zero_inds);


%%% Interpolate the Grid
disp('Writing Data')


    slopegrid_x = opt_x;
    slopegrid_y = opt_y;
    slopegrid = opt_angle;

disp(['Line Complete'])


end