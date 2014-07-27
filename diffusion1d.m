% Parameters
% ==========
%    kappa : the diffusivity coefficient
%    x_int : the space range [a, b]
%    t_int : the time interval [t0, tfinal] 
%
%    u_init : a function giving the initial temperature from points a to b.
%    u_bndry : a function giving the two boundary temperatures, at time t
%
%    n_x : the number of points into which we will divide [a, b]
%    n_t : the number of points into which we will divide [t0, tfinal]
%
% Temporary Values
% =============
%   dt: Delta time, its a small increment in time. It is calculated by
%   diving the t_int (time interval) into n-1 divisions
%   h:  A small increment in the position, calculated by diving the x_int
%   (x interval) into n-1 divisions
%   correct_nt: The correct n_t value that the user is prompted if the
%   ratio of kappa*dt/h^2 is greater than 0.5 which would cause the system
%   to diverge.
%   u_bndry_values: Contains the boundary temperatures at the boundary x
%   points.
%   k: its the counter, used to iterate through the columns of the u
%   matrix.
% 
% Return Values
% =============
%    x : returns a vector that contains nx points from a to b.
%    t : returns a vector that contains nt points from tinitial .
%    u : returns a matrix containing the temperature at a given time and a x point.

function [x, t, u] = diffusion1d( kappa, x_int, nx, t_int, n_t, u_init, u_bndry )
       
    % Error Checking
    % ==============
    % check if the ratio of kappa*dt/h^2 is greater than 0.5, then tell 
    % the user to divide the time range further, by giving a value for n_t
    % that would satisfy the condition of the ratio of kappa*dt/h^2  being 
    % less than 0.5.

    dt = (t_int(end) - t_int(1))/ (n_t-1);
    h = ( x_int(end) - x_int(1) )/(nx -1);
    correct_nt = floor((2*kappa*( t_int(end) - t_int(1) )/ h^2) + 2);
    if ((kappa*(dt))/ (h^2) >= 0.5)
        throw ( MException( 'MATLAB:invalid_argument', ...
            'The ratio kappa*dt/h^2 = %f >= 0.5, consider using n_t = %d',(kappa*(dt))/ (h^2),...
            correct_nt));
    end

    % Initialization
    % ==============
    % Create a nx by n_t matrix containing the boundary values and the 
    % initial temperature values.
    
    x = linspace(x_int(1), x_int(end), nx);
    t = linspace(t_int(1), t_int(end), n_t);
    u = zeros (nx, n_t);
    u_bndry_values = u_bndry(linspace(t_int(1),t_int(end),n_t-1));
    u(1,2:end) = u_bndry_values (1, 1:end);
    u(end,2:end) = u_bndry_values (2, 1:end);
    u(:,1) = u_init(linspace(x_int(1),x_int(end),nx));
    
    % Solving
    % =======
    % Using the second difference from the previous column and applying the
    % formula from the slides the values in the current column are calculated .
    for k = 1: n_t-1
      u(2:end-1,k+1) = u(2:end-1, k) + (kappa*dt/h^2)*(diff(u(:,k),2));
    end
   
       
    
end
    
