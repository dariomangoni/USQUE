% USQUE filter

clear all
clc

load 'DataTest.mat'
% these data files are loaded from DataTest
acc_vera = acc_sensore_noisy;
% adding noise to sensors signal
GPS_meas = GPS_meas + 0.000001*randn(size(GPS_meas));
sonic_altitude = sonic_altitude + 0.05*randn(size(sonic_altitude));
acc_sensore_noisy = acc_sensore_noisy + 0.01*randn(size(acc_sensore_noisy))+0.025;
omega_sensore_noisy = omega_sensore_noisy + 0.001*randn(size(acc_sensore_noisy))+0.005;
global gps_init;

stp_fin = 3000;

%% Initalize filter
x_posteriori=zeros(16,stp_fin);
x_posteriori(:,1)=[1;0;0;0;... % initial quaternion
                   0;0;0;...   % initial omega bias
                   0;0;0;...   % initial position
                   0;0;0;...   % initial speed
                   0;0;0];     % initial acceleration bias

Px_pos = 0.4; % initial variance of state
Px_vel = 0.1;
Px_acc = 0.00001;

Rv_pos = 1e-2;
Rv_vel = 1e-8;
Rv_acc = 1e-6;

Rn_pos = 1;
Rn_vel = 10;

Px = diag([0.25,0.25,0.25,   0.25,0.25,0.25,   Px_pos,Px_pos,Px_pos   Px_vel,Px_vel,Px_vel,   Px_acc,Px_acc,Px_acc]);
Rv = diag([1e-12,1e-12,1e-12,   1e-12,1e-12,1e-12,   Rv_pos,Rv_pos,Rv_pos,   Rv_vel,Rv_vel,Rv_vel,   Rv_acc,Rv_acc,Rv_acc]);
Rn = diag([ 0.0001,0.0001,0.0001,   Rn_pos,Rn_pos,Rn_pos,  Rn_vel,Rn_vel]);


%% Filtering
tic
for stp=2:stp_fin
    [x_posteriori(:,stp), Px] = filtroCrassidis ( x_posteriori(:,stp-1), Px, Rv, Rn,  [omega_sensore_noisy(stp,:)'; acc_sensore_noisy(stp,:)'; vettore_magn; GPS_meas(stp,:)'; GPS_vel(stp,:)'; sonic_altitude(stp)] ,dt);
end
toc

close all
clear Px Rv Rn

coord_percorso = zeros(3,stp_fin);
for stp=1:stp_fin
    coord_percorso(:,stp) = GPS2coord(gps_init,GPS_meas(stp,:));
end

%pos
figure('Name','Cinematica','NumberTitle','Off','OuterPosition',[80 50 1480 850]);
subplot(3,3,1); title('pos X'); hold on;
plot(coord_percorso(1,:),'r');
plot(x_posteriori(8,:),'b','MarkerSize',0.5);
plot(pos(1:stp_fin),'g');
legend('sensor','filtered','real');

subplot(3,3,2); title('pos Y'); hold on;
plot(coord_percorso(2,:),'r');
plot(x_posteriori(9,:),'b','MarkerSize',0.5);
plot(zeros(stp_fin,1),'g');
legend('sensor','filtered','real');

subplot(3,3,3); title('pos Z'); hold on;
plot(sonic_altitude(1:stp_fin),'r');
plot(x_posteriori(10,:),'b','MarkerSize',0.5);
plot(zeros(stp_fin,1),'g');
legend('sensor','filtered','real');

% acc 
subplot(3,3,7); title('acc X'); hold on;
plot(acc_sensore_noisy(1:stp_fin,1),'b');
plot(acc_sensore_noisy(1:stp_fin,1)'-x_posteriori(14,1:stp_fin),'g');
plot(acc_vera(1:stp_fin,1),'r');
legend('sensor','filtered','real');

subplot(3,3,8); title('acc Y'); hold on;
plot(acc_sensore_noisy(1:stp_fin,2),'b');
plot(acc_sensore_noisy(1:stp_fin,2)'-x_posteriori(15,:),'g');
plot(acc_vera(1:stp_fin,2),'r');
legend('sensor','filtered','real');

subplot(3,3,9); title('acc Z'); hold on;
plot(acc_sensore_noisy(1:stp_fin,3),'b');
plot(acc_sensore_noisy(1:stp_fin,3)'-x_posteriori(16,1:stp_fin),'g');
plot(acc_vera(1:stp_fin,3),'r');
legend('sensor','filtered','real');


subplot(3,3,4); title('vel X'); hold on;
plot(vel(1:stp_fin,1),'r');
plot(x_posteriori(11,1:stp_fin),'c');
legend('real','filtered');

subplot(3,3,5); title('vel Y'); hold on;
plot(vel(1:stp_fin,2),'r');
plot(x_posteriori(12,:),'c');
legend('real','filtered');


subplot(3,3,6); title('vel Z'); hold on;
plot(vel(1:stp_fin,3),'r');
plot(x_posteriori(13,1:stp_fin),'c');
legend('real','filtered');
