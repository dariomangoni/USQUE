
function [state_new, Px] = filtroCrassidis ( state_old, Px, Q, R, measurement, dt )

global gps_init;
state_new=state_old;

a = 1;
f = 2*(a+1);
n = 15;
lambda = 1;
n_obs = 8;

%% Allocazione matrici
q      = zeros( 4, 2*n+1 );
sigmaY = zeros( n_obs, 2*n+1 );
dq     = zeros( 4, 2*n+1 );

%% Inizializzazione matrici
omega_bias = state_old(5:7);
pos = state_old(8:10);
vel = state_old(11:13);
acc_bias = state_old(14:16);


omega_meas = measurement(1:3);
acc_meas = measurement(4:6);
magn_meas = measurement(7:9);
gps_meas = measurement(10:12);
gps_vel = measurement(13:14);
gps_diff = GPS2coord(gps_init, gps_meas);
gps_vel = GPS2vel(gps_vel);
sonic_altitude = measurement(15);

%% %%%%%%%%%%%%%%%%%%%%% Predict part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calcolo varianza processo
PQ = chol((n+lambda)*(Px+Q))';

%% creo il sigmaset (5)
sigmaX(:,1) = [0;0;0;omega_bias;pos;vel;acc_bias];
for s=2:n+1
    sigmaX(:, s     ) = sigmaX(:,1) + PQ(:,s-1);
    sigmaX(:, s + n ) = sigmaX(:,1) - PQ(:,s-1);
end

%% creo l'error-quaternion dq (33)
dq(:,1) = [1;0;0;0];
for s=2:(2*n+1)
    normsigmaX2 = sigmaX(1:3,s)'*sigmaX(1:3,s);
    dq(1,s) = (-a*normsigmaX2 + f*sqrt( f^2 + (1-a^2)*normsigmaX2 ) ) / ( f^2 + normsigmaX2 );
    dq(2:4,s) = (a + dq(1,s)) * sigmaX(1:3,s) / f;
end

%% creo sigmapoint quaternions che verranno propagati(32)
q(:,1) = state_old(1:4);
for s=2:(2*n+1)
    q(:,s) = quaternione_moltiplicazione( dq(:,s), q(:,1) );
end

%% propagated quaternions and state (34) - process model
for s=1:(2*n+1)
%     q(:,s) = quaternione_aggiornamento( q(:,s), ( omega_meas - sigmaX(4:6,s) ), dt );
%     q(:,s) = quaternione_aggiornamento( q(:,s), ( omega_meas  ), dt );
      [q(:,s), sigmaX(:,s)] = process_model( q(:,s), sigmaX(:,s), dt, acc_meas, omega_meas);
end

%% propagated error quaternions (36) - scarti dal valore medio
for s= 1:(2*n+1)
    dq(:,s) = quaternione_moltiplicazione( q(:,s) , q(:,1).*[1;-1;-1;-1] );
end

%% propagated sigma-points (37,38)
sigmaX(1:3,1) = [0;0;0];
for s=2:(2*n+1)
    sigmaX(1:3,s) = f * dq(2:4,s)/( a + dq(1,s) );
end

%% predicted mean (7) - può essere messo subito dopo a (34)
x_priori = ( lambda * sigmaX(:,1) + 0.5 * sum( sigmaX(:,2:end) ,2) ) / ( n + lambda );

%% predicted covariance (8)
sigmaXscarto(:,1) = sigmaX(:,1) - x_priori;
Px_priori = lambda * sigmaXscarto(:,1) * sigmaXscarto(:,1)';
for s = 2:(2*n+1)
    sigmaXscarto(:,s) = sigmaX(:,s) - x_priori;
    Px_priori = Px_priori + 0.5 * sigmaXscarto(:,s) * sigmaXscarto(:,s)';
end
Px_priori = Px_priori / ( n + lambda ) + Q;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update part %%%%%%%%%%%%%%%%%%%%%%%%%%

%% observation model (10+43)
for s = 1:(2*n+1)
%     sigmaY(:,s) = earth_to_body( [0;0;-9.8], q(:,s) ); % vedi (43)
    sigmaY(:,s) = observation_model( q(:,s), sigmaX(:,s) );
end

%% observation mean (9)
y_priori = ( lambda * sigmaY(:,1) + 0.5 * sum( sigmaY(:,2:end), 2 ) )/ ( n+lambda );

%% meas e output covariance (11,12)
sigmaYscarto(:,1) = sigmaY(:,1) - y_priori;
Py = lambda * sigmaYscarto(:,1) * sigmaYscarto(:,1)';
for s = 2:(2*n+1)
    sigmaYscarto(:,s) = sigmaY(:,s) - y_priori;
    Py = Py + 0.5 * sigmaYscarto(:,s) * sigmaYscarto(:,s)';
end
Py = Py / ( n + lambda );
Pv = Py + R;

%% cross-correlation matrix (13)
Pxy = lambda * sigmaXscarto(:,1) * sigmaYscarto(:,1)';
for s =2:(2*n+1)
    Pxy =  Pxy + 0.5 * sigmaXscarto(:,s) * sigmaYscarto(:,s)';
end
Pxy = Pxy / ( n+lambda );

%% updating state vector
K = Pxy/Pv;
innov = [magn_meas;gps_diff(1:2);sonic_altitude;gps_vel] - y_priori;
x_posteriori = x_priori + K*innov;

normdp2 = x_posteriori(1:3)'*x_posteriori(1:3);
dq(1,1) = ( -a * normdp2 + f * sqrt( f^2 + (1-a^2) * normdp2 ) ) / ( f^2 + normdp2 );
dq(2:4,1) = x_posteriori(1:3) * ( a + dq(1,1) ) / f;

state_new(1:4) = quaternione_moltiplicazione( dq(:,1) , q(:,1) );
state_new(5:16) = x_posteriori(4:15);

Px = Px_priori - K*Pv*K';

end
