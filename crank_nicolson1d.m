function [x, t, u] = crank_nicolson1d( kappa, x_int, nx, t_int, n_t, u_init, u_bndry )
 
  
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
%
% 
% Return Values
% =============
%    x : returns a vector that contains nx points from a to b.
%    t : returns a vector that contains nt points from tinitial .
%    u : returns a matrix containing the temperature at a given time and a x point.

% Error Checking
% ==============

    if ~all( size( x_int ) == [2, 1] )  % ensures that there are 2 boundary x values in the form of a 2 x 1 column vector
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument x_int is not a 2-dimensional column vector' ) );
    end
    
    if ( nx < 4) %ensures that there are enough nx points to calculate the values at the insulated boundary condition
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument nx is smaller than 4' ) );
    end
    
    if ~all( size( t_int ) == [2, 1] )  % ensures that there are 2 boundary u values in the form of a 2 x 1 column vector
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_int is not a 2-dimensional column vector' ) );
    end

% Initialization
% ==============
% Create a nx by n_t matrix containing the boundary values and the 
% initial temperature values.
      
    dt = (t_int(end) - t_int(1))/ (n_t-1);
    h = ( x_int(end) - x_int(1) )/(nx -1);
    if ((kappa*(dt))/ (h^2) >= 0.5)
        warning( 'MATLAB:questionable_argument', ...
        'the kappa value is sub-optimal, there may be decaying oscillations')
    end
 
    x = linspace(x_int(1), x_int(end), nx);
    t = linspace(t_int(1), t_int(end), n_t);
    u = zeros (nx, n_t);
    u_bndry_values = u_bndry(t);
  
    u(1,2:end) = u_bndry_values (1, 2:end);
    u(end,2:end) = u_bndry_values (2, 2:end);
    u(:,1) = u_init(x);
    r = kappa*dt/h^2;
    
    %Calculate M1 , M2 and M3 which are later used to calculate M.
    %M1 , M2 and M3 are diagonal matrices
    M1 = diag(2*(r+1)*ones(nx-2,1));
    M2 = diag((-r)*ones(nx-3, 1),1);
    M3 = diag((-r)*ones(nx-3, 1),-1);
      
    for k = 1:n_t-1
        M = M1 + M2 + M3;
        % Calculate b vector using the formula in the presentation
        b = 2*u(2:end-1, k) + r*(diff(u(:,k),2));
    
        % If the top boundary is an insulator the 
        % following calculations are performed on matrix M and b
        if(isnan(u(1,k+1))) 
            M(1,1) = 2 + 2*r/3;
            M(1,2) = -2*r/3;
        else %the first value of b is changed according to the formula
            b(1) = b(1) +  r*u_bndry_values(1,k+1);
        end
        % If the lower boundary is an insulator the 
        % following calculations are performed on matrix M and b
        if(isnan(u(end,k+1)))
            M(end,end) = 2 + 2*r/3;
            M(end, end-1) = -2*r/3;
        else% the last value of b is changed according to the formula
            b(end,1) = b(end,1) + r*u_bndry_values(2,k+1);
        end
        % The matrix solution is calulated by dividing b by M
        % and is set to the apporiate matrix
        u_answer = M\b;
        u(2:end-1,k+1) = u_answer;
        
        if( isnan(u(1,k+1)) )
            u(1,k+1) = 4*u(2,k+1)/3 -u(3,k+1)/3;
        end
        
        if( isnan(u(end,k+1)) )
            u(end,k+1) = 4*u(end-1,k+1)/3 -u(end-2,k+1)/3;
        end
        
  
    end
       
end
