function [vel] = GPS2vel (gps_vel)
vel_abs = gps_vel(1)/3.6; % [m s^-1] speed modulus
vel_ang = gps_vel(2)*pi/180; % speed direction in radians

vel = zeros(2,1); %initialize vector
vel(1) = vel_abs * sin(vel_ang); % speed along global x axis
vel(2) = vel_abs * cos(vel_ang); % speed along global y axis
end