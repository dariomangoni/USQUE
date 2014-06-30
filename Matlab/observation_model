function [measurement_attended] = observation_model ( quat, sigmaX )

vettore_magn = [0.48528001;0.018560000;-0.87415999]; % magnetic vector measured when quadcopter points to north
measurement_attended = zeros(8,1); % initialize vector

measurement_attended(1:3) = earth_to_body( vettore_magn , quat ); % magnetometer measurement attended
measurement_attended(4:6) = sigmaX(7:9); % GPS position measurement attended
measurement_attended(7:8) = sigmaX(10:11); % GPS speed measurement attended
end
