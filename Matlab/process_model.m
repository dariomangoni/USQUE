function [quat, sigmaX] = process_model (quat, sigmaX, dt, acc_meas, omega_meas)

acc_vero = acc_meas - sigmaX(13:15); % subtract acceleration bias from acceleration
omega_vero = omega_meas - sigmaX(4:6); % subtract oemga bias from omega

%calculate rotation matrix from quaternion
rotation_matrix=2*[0.5-quat(3)^2-quat(4)^2,         quat(2)*quat(3)-quat(1)*quat(4), quat(2)*quat(4)+quat(1)*quat(3);...
                   quat(2)*quat(3)+quat(1)*quat(4), 0.5-quat(2)^2-quat(4)^2,         quat(3)*quat(4)-quat(1)*quat(2);...
                   quat(2)*quat(4)-quat(1)*quat(3), quat(3)*quat(4)+quat(1)*quat(2), 0.5-quat(2)^2-quat(3)^2        ];

vel_1 = rotation_matrix * acc_vero + [0;0;-9.81]; % acceleration in global reference

%% Update variables
sigmaX(10:12) = sigmaX(10:12) + vel_1*dt; % speed update
sigmaX(7:9) = sigmaX(7:9) + sigmaX(10:12)*dt; % position update
quat = quaternione_aggiornamento( quat, omega_vero, dt ); % quaternion update
end
