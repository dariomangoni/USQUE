function [delta_pos] = GPS2coord ( posINIT, posNEW )
% from coordinates to position (in meters) from a initial point
delta_pos=zeros(3,1);
earth_radius = 6362697;
delta_pos(1) = ( posNEW(2) - posINIT(2) )*pi/180*earth_radius*cosd(posNEW(2));
delta_pos(2) = ( posNEW(1) - posINIT(1) )*pi/180*earth_radius;
delta_pos(3) = posNEW(3) - posINIT(3);
end
