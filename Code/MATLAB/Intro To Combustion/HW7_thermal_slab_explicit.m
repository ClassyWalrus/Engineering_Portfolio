%Solves the 2D heat equation with an explicit finite difference scheme
%dT/dt = k * d^2T/dt^2
%%
clear

%Physical parameters
k=2; % heat constant
L=10; %length of tube
T_initial = 300;    % initial temperature in K
T_boundary = 500;   % boundary temperature in K





% discretization
N=10; % number of elements in space
Nt=20;
x_vec=linspace(0,L,N);
dx=x_vec(2)-x_vec(1);
alpha = [0.5 , 1, 2];
dt= alpha*(dx^2)/(2*k);
t_vec=0:dt:Nt;

%initialization
T_mat=zeros(length(x_vec),length(t_vec));
T_imp=zeros(length(x_vec),length(t_vec));
T_mat(:,1)=300;
T = T_initial * ones(N, 1); % initial condition
T(1) = T_boundary;  % left boundary condition
T(end) = T_boundary; % right boundary condition




%all zero

%boundary conditions
T_mat(1,:)=500;
T_mat(end,:)=500;
for j = 1:length(dt)

    A = zeros(N, N);
    for i = 2:N-1
        A(i,i-1) = -alpha(j); % lower diagonal
        A(i,i) = 1 + 2 * alpha(j); % main diagonal
        A(i,i+1) = -alpha(j); % upper diagonal
    end
    A(1,1) = 1;        % Boundary condition at the left wall
    A(end,end) = 1;    % Boundary condition at the right wall


    %solving
    for tdx=1:length(t_vec)-1
        for idx=2:length(x_vec)-1   
            T_mat(idx,tdx+1)=T_mat(idx,tdx)+k*dt(j)/(dx^2)*(T_mat(idx+1,tdx)-2*T_mat(idx,tdx)+T_mat(idx-1,tdx));
        end

% Right-hand side vector
        b = T;
        b(2:N-1) = T(2:N-1) + alpha(j) * (T(1:N-2) - 2*T(2:N-1) + T(3:N));

        % Solve the linear system A*T_new = b
        T(j) = A \ b;
    end
    



    Ex_C = T_mat();
    Imp_C = T_imp();
    
    
    
    %plotting
    figure(1)
    [tt,xx]=meshgrid(t_vec,x_vec);
    mesh(xx,tt,T_mat)
    xlabel('x coordinate (m)')
    ylabel('Time(s)')
    zlabel('Temperature (K)')
    
    
    % figure(2)
    % mesh(xx,tt,T_imp)
    
 end
    %%
    
%     % Initialize temperature array
%     T = T_initial * ones(N, 1); % initial condition
%     T(1) = T_boundary;  % left boundary condition
%     T(end) = T_boundary; % right boundary condition
% 
%     % Coefficient matrix for implicit method
%     A = zeros(N, N);
%     for i = 2:N-1
%         A(i,i-1) = -alpha; % lower diagonal
%         A(i,i) = 1 + 2 * alpha; % main diagonal
%         A(i,i+1) = -alpha; % upper diagonal
%     end
%     A(1,1) = 1;        % Boundary condition at the left wall
%     A(end,end) = 1;    % Boundary condition at the right wall
% 
%     % Matrix to store temperature data at each time step for plotting
%     T_matrix = zeros(N, Nt); % Each column corresponds to temperature at a time step
% 
%     % Time-stepping loop
%     for n = 1:Nt
%         % Right-hand side vector
%         b = T;
%         b(2:N-1) = T(2:N-1) + alpha * (T(1:N-2) - 2*T(2:N-1) + T(3:N));
% 
%         % Solve the linear system A*T_new = b
%         T_new = A \ b;
% 
%         % Update the temperature array
%         T = T_new;
% 
%         % Store the temperature profile at the current time step
%         T_matrix(:, n) = T;
%     end

% % Create meshgrid for plotting
% figure(2)
% x = linspace(0, L, N);        % Spatial grid (0 to L)
% t = linspace(0, Nt*dt, Nt);    % Time grid (0 to final time)
% mesh(x,t,T_matrix')
% xlabel('x coordinate (m)')
% ylabel('Time(s)')
% zlabel('Temperature (K)')
% title('Temperature Distribution Over Time (3D Meshgrid)');