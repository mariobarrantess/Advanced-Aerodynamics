%% AERODYNAMICS 2.- PROJECT 2 2025-2026 - WING #19
% MARIO BARRANTES MELLADO 100471762
% GUILLERMO MARÍN COBO 100475752

clc; clear all; close all;
%% Parameters of Wing & ISA

% Parameters of the wing
Sw = 55;             % Surface of the wing [m^2] - from statement
AR = 6.14;           % Area Ratio = b^2/S - from XFLR5
b = sqrt(AR*Sw);     % Wingspan [m]

% Ni = Number of panels in each direction
Ny = 80;
Nx = 40;

% ISA PARAMETERS
M = 22000;             % Mass of the aircraft [kg] - from statement
g = 9.81;              % Acceleration due to gravity force [m/s^2]
CD0 = 0.01;            % Parasitic Drag coefficient - from statement
h = 11000;             % Flying altitude [m] - from statement
T0 = 288.15;           % Ambient temperature @ h=0 [K]
A = 0.0065;            % Lapse rate - ΔT/m [K/m] 
R = 287.05;            % Ideal gas constant [J·kg/K]
T = T0 - A * h;        % Formula to compute temperature @ certain altitude [K]
P = 101325 * (1-(A*h)/T0)^(g/(A * R)); % Formula to compute pressure @ certain altitude [Pa]
rho = P/(R*T);         % Density of air @ h=11000m [kg/m^3] - rho = P/(R*T)
e = 1;                 % Oswald Efficiency Factor
gamma_isa = 1.4 ;          % Ideal gas ratio
W = M * g;             % Weight of aircraft [N]
a = sqrt(gamma_isa*R*T);   % Speed of sound [m/s]

b_0 = b;   % Keep actual value of span [m]

%% WING GEOMETRY for EIW - M_inf=[0,0.1,0.3,0.5,0.7]
M_inf = [0,0.1,0.3,0.5,0.7];

beta = zeros(1, length(M_inf));
b_eiw = zeros(1, length(M_inf));
b_kink_eiw = zeros(1, length(M_inf));
Sw_eiw =  zeros(1, length(M_inf));

colors = lines(length(M_inf)); % At the start, outside the loop
handles = zeros(1, length(M_inf)); % Preallocate for legend

for k=1:length(M_inf)

beta(k) = sqrt(1-(M_inf(k)^2));  % Prandtl-Glauert Correction Factor
b_eiw (k) =  b_0*beta(k);          % Span of Equivalent Incompressible Wing
Sw_eiw (k) = Sw*beta(k);           % Surface of Equivalent Incompressible Wing

%hay que repetir todo el proceso del project1 para el equivalent incompressible wing

% only compress y-coordinates
b = b_eiw(k);             % Wingspan [m]
c_root = b_0/3.5;        % Chord @ root [m]
% b_kink = c_root * 1.6 ;   % Span @ kink [m]

b_kink = (1.6/3.5)*b; % Span @ kink - f(b) not f(c_root)
c_tip = c_root * 0.2;    % Chord @ tip [m]
lambda = 28;             % Swept angle [deg]
lambda_eiw = atan(tand(lambda)/beta(k)) * (180/pi); % Equivalent sweep angle in degrees

% Coordinates of characteristic points of the wing
tip = [0,0] ;
b_kink_top = [(b_kink/2) * tand(lambda_eiw) , b_kink/2];
b_kink_bottom = [c_root , b_kink/2];
c_tip_top = [(b/2)*tand(lambda_eiw) , b/2];
c_tip_bottom = [(b/2)*tand(lambda_eiw) + c_tip , b/2];
c_root_profile = [c_root , 0];

% Define the 5 points of the wing outline - 3 sections where geometry changes: Tip, Kink & Root
x_wing_profile = [tip(1), b_kink_top(1), c_tip_top(1), c_tip_bottom(1), b_kink_bottom(1), c_root_profile(1)]; %, tip(1)
y_wing_profile = [tip(2), b_kink_top(2), c_tip_top(2), c_tip_bottom(2), b_kink_bottom(2), c_root_profile(2)]; %, tip(2)
y_wing_profile_left = -y_wing_profile;

%%% Wing geometry and panel generation%%
% Storing panels points in a matrix, the same matrix for all the wing
PanelPoints = zeros((Ny/2)+1, Nx+1, 2); 
PanelPoints_full = zeros(Ny+ 1, Nx+1, 2) ;

% Defining the leading edge coordinates
Leading_edge_root_kink = [tip(1), tip(2); b_kink_top(1), b_kink_top(2)];         % LE Root-Kink
Leading_edge_kink_tip = [b_kink_top(1), b_kink_top(2); c_tip_top(1), c_tip_top(2)]; % LE Kink-Tip

% Defining the trailing edge coordinates
Trailing_edge_root_kink = [c_root_profile(1), c_root_profile(2); b_kink_bottom(1), b_kink_bottom(2)]; % TE Root-Kink
Trailing_edge_kink_tip = [b_kink_bottom(1), b_kink_bottom(2); c_tip_bottom(1), c_tip_bottom(2)];   % TE Kink-Tip

dy = b_eiw(k)/Ny;          % Width of each panel
dx = linspace(0, 1, Nx+1); % Length of each panel

for j = 1:(Ny/2) + 1 
    y_j = (j-1) * dy; 

    if y_j > b/2
        y_j = b/2;
    end
    % local_x_le = 0.0;
    % local_x_te = 0.0;

    % Interpolate LE and TE x-coordinates depending if we are before or after the kink
    if y_j <= b_kink_top(2) % Root-kink
        local_x_le = interp1(Leading_edge_root_kink(:,2), Leading_edge_root_kink(:,1), y_j);  % LE x-coordinate root-kink
        local_x_te = interp1(Trailing_edge_root_kink(:,2), Trailing_edge_root_kink(:,1), y_j); % TE x-coordinate
    else % Kink-tip
        local_x_le = interp1(Leading_edge_kink_tip(:,2), Leading_edge_kink_tip(:,1), y_j);  % LE x-coordinate kink-tip
        local_x_te = interp1(Trailing_edge_kink_tip(:,2), Trailing_edge_kink_tip(:,1), y_j); % TE x-coordinate kink-tip
    end

    % Calculates the chord length of each panel
    local_chord = local_x_te - local_x_le;

    % Creates all chordwise panel vertices a this y-coordinate
    for i = 1:(Nx+1)
        PanelPoints(j, i, 1) = local_x_le + dx(i) * local_chord; % x-coordinate of panel i,j
        PanelPoints(j, i, 2) = y_j;         % y-coordinate of panel i,j
    end
end

% Mirroring half of the wing to get the complete wing plot
PanelPoints_full(1:Ny/2+1, :, 2) = -PanelPoints(end:-1:1,:,2) ;
PanelPoints_full(Ny/2+2 : Ny + 1, :, 2) = PanelPoints(2:end,:,2);
PanelPoints_full(1:Ny/2+1, :, 1) = PanelPoints(end:-1:1,:,1);
PanelPoints_full(Ny/2+2 : Ny + 1, :, 1) = PanelPoints(2:end,:,1) ;


% figure; % Plotting the wing surface

col = colors(k, :); % Pick a color for this Mach number

h = plot(y_wing_profile, -x_wing_profile, 'Color', col, 'LineWidth', 1.5);
hold on
plot(y_wing_profile_left, -x_wing_profile, 'Color', col, 'LineWidth', 1.5);
hold on

handles(k) = h; % Store handle for legend

% % Plotting panels division inside the wing
% for j = 1:Ny+1
    % plot(PanelPoints_full(j, :, 2), -PanelPoints_full(j, :, 1), 'b', 'LineWidth', 0.5);
    % hold on
% end
% for i = 1:Nx+1
    % plot(PanelPoints_full(:, i, 2), -PanelPoints_full(:, i, 1), 'b', 'LineWidth', 0.5);
    % hold on
% end
axis equal % Centering axis
hold on
end

title('Equivalent Incompressible Wings for $M_{\infty}=[0, 0.1,0.3,0.5,0.7]$', Interpreter='latex', FontWeight='bold', FontSize=14)
legendStrings = arrayfun(@(m) sprintf('M{\\infty} = %.1f', m), M_inf, 'UniformOutput', false);
legend(handles, legendStrings, Interpreter= 'latex');
xlabel('y')
ylabel('x')
axis equal
grid on
yticks = get(gca, 'YTick');          % This grabs MATLAB's current y-axis ticks
yticklabels = -yticks;               % Make all tick labels positive
set(gca, 'YTickLabel', yticklabels); % This replaces label -3 with 3, etc.
hold off

%% PANEL DIVISION FOR EIW - M_inf=[0,0.1,0.3,0.5,0.7]
%k = 1; % To plot panel division for incompressible wing (M_inf = 0) - same as in project 1
% b = sqrt(AR*Sw);     % Wingspan [m]

for k = 1:length(M_inf)

% para probar si funciona bien, ploteo solo Mach menor (M=0) y mayor (M=0.7)
% if k==1 || k==5 % only entering loop (plotting wing and panel divisions + control points) for M=0 and M=0.7

beta(k) = sqrt(1-(M_inf(k)^2));  % Prandtl-Glauert Correction Factor
b_eiw (k) =  b_0*beta(k);          % Span of Equivalent Incompressible Wing

%hay que repetir todo el proceso del project1 para el equivalent incompressible wing

% only compress y-coordinates
b = b_eiw(k);             % Wingspan [m]
c_root = b_0/3.5;        % Chord @ root [m]
% b_kink = c_root * 1.6 ;   % Span @ kink [m]
b_kink = (1.6/3.5)*b; % Span @ kink - f(b) not f(c_root)
c_tip = c_root * 0.2;    % Chord @ tip [m]
lambda = 28;             % Swept angle [deg]
lambda_eiw = atan(tand(lambda)/beta(k)) * (180/pi); % Equivalent sweep angle in degrees


% Coordinates of characteristic points of the wing
tip = [0,0] ;
b_kink_top = [(b_kink/2) * tand(lambda_eiw) , b_kink/2];
b_kink_bottom = [c_root , b_kink/2];
c_tip_top = [(b/2)*tand(lambda_eiw) , b/2];
c_tip_bottom = [(b/2)*tand(lambda_eiw) + c_tip , b/2];
c_root_profile = [c_root , 0];

% Define the 6 points of the wing outline - 3 sections where geometry changes: Tip, Kink & Root
x_wing_profile = [tip(1), b_kink_top(1), c_tip_top(1), c_tip_bottom(1), b_kink_bottom(1), c_root_profile(1)]; %, tip(1)
y_wing_profile = [tip(2), b_kink_top(2), c_tip_top(2), c_tip_bottom(2), b_kink_bottom(2), c_root_profile(2)]; %, tip(2)
y_wing_profile_left = -y_wing_profile;


%%% Wing geometry and panel generation%%
% Storing panels points in a matrix, the same matrix for all the wing
PanelPoints = zeros((Ny/2)+1, Nx+1, 2); 
PanelPoints_full = zeros(Ny+ 1, Nx+1, 2) ;

% Defining the leading edge coordinates
Leading_edge_root_kink = [tip(1), tip(2); b_kink_top(1), b_kink_top(2)];         % LE Root-Kink
Leading_edge_kink_tip = [b_kink_top(1), b_kink_top(2); c_tip_top(1), c_tip_top(2)]; % LE Kink-Tip

% Defining the trailing edge coordinates
Trailing_edge_root_kink = [c_root_profile(1), c_root_profile(2); b_kink_bottom(1), b_kink_bottom(2)]; % TE Root-Kink
Trailing_edge_kink_tip = [b_kink_bottom(1), b_kink_bottom(2); c_tip_bottom(1), c_tip_bottom(2)];   % TE Kink-Tip

dy = b_eiw(k)/Ny;               % Width of each panel
dx = linspace(0, 1, Nx+1); % Length of each panel

for j = 1:(Ny/2) + 1 
    y_j = (j-1) * dy; 

    if y_j > b/2
        y_j = b/2;
    end
    % local_x_le = 0.0;
    % local_x_te = 0.0;

    % Interpolate LE and TE x-coordinates depending if we are before or after the kink
    if y_j <= b_kink_top(2) % Root-kink
        local_x_le = interp1(Leading_edge_root_kink(:,2), Leading_edge_root_kink(:,1), y_j);  % LE x-coordinate root-kink
        local_x_te = interp1(Trailing_edge_root_kink(:,2), Trailing_edge_root_kink(:,1), y_j); % TE x-coordinate
    else % Kink-tip
        local_x_le = interp1(Leading_edge_kink_tip(:,2), Leading_edge_kink_tip(:,1), y_j);  % LE x-coordinate kink-tip
        local_x_te = interp1(Trailing_edge_kink_tip(:,2), Trailing_edge_kink_tip(:,1), y_j); % TE x-coordinate kink-tip
    end

    % Calculates the chord length of each panel
    local_chord = local_x_te - local_x_le;

    % Creates all chordwise panel vertices a this y-coordinate
    for i = 1:(Nx+1)
        PanelPoints(j, i, 1) = local_x_le + dx(i) * local_chord; % x-coordinate of panel i,j
        PanelPoints(j, i, 2) = y_j;         % y-coordinate of panel i,j
    end
end

% Mirroring half of the wing to get the complete wing plot
PanelPoints_full(1:Ny/2+1, :, 2) = -PanelPoints(end:-1:1,:,2) ;
PanelPoints_full(Ny/2+2 : Ny + 1, :, 2) = PanelPoints(2:end,:,2);
PanelPoints_full(1:Ny/2+1, :, 1) = PanelPoints(end:-1:1,:,1);
PanelPoints_full(Ny/2+2 : Ny + 1, :, 1) = PanelPoints(2:end,:,1) ;


figure; % Plotting the wing surface

col = colors(k, :); % Pick a color for this Mach number

h = plot(y_wing_profile, -x_wing_profile, 'Color', col, 'LineWidth', 1.5);
hold on
plot(y_wing_profile_left, -x_wing_profile, 'Color', col, 'LineWidth', 1.5);
hold on

handles(k) = h; % Store handle for legend

% Plotting panels division inside the wing
for j = 1:Ny+1
    plot(PanelPoints_full(j, :, 2), -PanelPoints_full(j, :, 1), 'b', 'LineWidth', 0.5);
    hold on
end
for i = 1:Nx+1
    plot(PanelPoints_full(:, i, 2), -PanelPoints_full(:, i, 1), 'b', 'LineWidth', 0.5);
    hold on
end
axis equal % Centering axis
hold on

%%% Computing and plotting Control Points
x_c = zeros(Nx, Ny); % Control point x-coordinates
y_c = zeros(Nx, Ny); % Control point y-coordinates

for i = 1:Nx 
    for j = 1:Ny 
        % x,y coordinates for the 4 corners of the (i,j) panel
        P1_coords = squeeze(PanelPoints_full(j,   i,   :)); % Inboard, LE
        P2_coords = squeeze(PanelPoints_full(j,   i+1, :)); % Inboard, TE
        P3_coords = squeeze(PanelPoints_full(j+1, i+1, :)); % Outboard, TE
        P4_coords = squeeze(PanelPoints_full(j+1, i,   :)); % Outboard, LE

        % control points @ exact center of the panel (same as in exercises of sheet 4)
        x_c(i,j) = (P1_coords(1) + P2_coords(1) + P3_coords(1) + P4_coords(1)) / 4;
        y_c(i,j) = (P1_coords(2) + P2_coords(2) + P3_coords(2) + P4_coords(2)) / 4;

        plot(y_c(i,j),-x_c(i,j),'r.', 'MarkerSize', 6) % Adding control points to plot of wing with panels division

    end
end

title(sprintf('Equivalent Incompressible Wing for $M_{\\infty} = %.2f$', M_inf(k)), Interpreter='latex', FontWeight='bold',FontSize=14)
xlabel('y (spanwise direction)', Interpreter= 'latex')
ylabel('x (chordwise direction)', Interpreter= 'latex')
axis on
grid off
limits = [min(y_wing_profile_left)-1,max(y_wing_profile)+1];
xlim (limits)
yticks = get(gca, 'YTick');      % This grabs MATLAB's current y-axis ticks
yticklabels = -yticks;       % Make all tick labels positive
set(gca, 'YTickLabel', yticklabels); % This replaces label -3 with 3, etc.
hold off

end

%% VLM FOR EIW TO GET dCL/dα for M_inf=[0,0.1,0.3,0.5,0.7] - from Project 1 - OK (revisado con Ignacio)
alpha_deg = [0, 5, 10, 15, 20] ; % Angle of attack [degrees]
alpha_rad = deg2rad(alpha_deg);  % Angle of attack [radians]

Na = length(alpha_rad);

CL_i = zeros(Na);
CL_c = zeros(Na);
L = zeros(Na);
CL_full_i = zeros(Na,Na);
CL_full_c = zeros(Na,Na);
dCL_dalpha_c_rad = zeros(1,Na);
dCL_dalpha_c_deg = zeros(1,Na);

colors = lines(length(M_inf));     % Vector for different colors for each Mach
handles = zeros(1, length(M_inf)); % Preallocate vector for legend

for i = 1:length(M_inf)

    [CL_i, ww_mat, Gamma, L, V_inf] = vlm (M_inf(i), PanelPoints_full, alpha_rad, Nx, Ny, x_c, y_c, b_eiw(i), rho, a, Sw_eiw(i));

    col2 = colors(i, :);      % Pick a color for plotting CL for this Mach number - same colors as EIWs

    CL_full_i(i,:) = CL_i;    % Store all values of CL - Each row is for each M_inf   -->  Equal for all M_inf
    CL_c = CL_i./beta(i);     % Undo PG contraction for CL
    CL_full_c(i,:) = CL_c;    % Store all values of CL - Each row is for each M_inf

    h2 = plot(alpha_deg, CL_c, 'Color', col2, 'LineWidth', 1.5);
    hold on
    handles(i) = h2;          % Store handles for legend

    dCL_dalpha_c_rad(i) = (CL_full_c(i, end) - CL_full_c(i, 1)) / (alpha_rad(end) - alpha_rad(1)); 
    dCL_dalpha_c_deg(i) = (CL_full_c(i, end) - CL_full_c(i, 1)) / (alpha_deg(end) - alpha_deg(1));

end

% Plot of CL vs alpha [deg]
title('$C_L$ vs $\alpha$ for $M_{\infty}=[0, 0.1,0.3,0.5,0.7]$', Interpreter='latex', FontWeight='bold', FontSize=14)
legendStrings = arrayfun(@(m) sprintf('M{\\infty} = %.1f', m), M_inf, 'UniformOutput', false);
legend(handles, legendStrings, Interpreter= 'latex');
xlabel('$\alpha$', Interpreter='latex', FontSize=14)
ylabel('$C_L$', Interpreter='latex', FontSize=14)
grid on
hold off

% Plot of dCL/dalpha [deg]
figure
plot (M_inf, dCL_dalpha_c_deg, 'o-', LineWidth=1.5)
title('$dC_L/d\alpha$ [deg] for $M_{\infty}=[0,0.1,0.3,0.5,0.7]$', Interpreter='latex', FontWeight='bold', FontSize=14)
xlabel('$M_{\infty}$', Interpreter='latex', FontSize=14)
ylabel('$dC_L/d\alpha$ [deg]', Interpreter='latex', FontSize=14)
grid on

% Plot of dCL/dalpha [rad]
figure
plot (M_inf, dCL_dalpha_c_rad, 'o-', Linewidth=1.5)
title('$dC_L/d\alpha$ [rad] for $M_{\infty}=[0,0.1,0.3,0.5,0.7]$', Interpreter='latex', FontWeight='bold', FontSize=14)
xlabel('$M_{\infty}$', Interpreter='latex', FontSize=14)
ylabel('$dC_L/d\alpha$ [rad]', Interpreter='latex', FontSize=14)
grid on



%%%%% COMPUTATION OF INDUCED DRAG FROM SAME PROCEDURE AS THE ONE USED FOR
%%%%% GETTING LIFT ON M_inf_D = [0,0.1,0.3,0.5,0.7] 
%% WING GEOMETRY for EIW - M_inf_D=[0.5-0.95]
M_inf_D = 0.5:0.05:0.95;           % Incoming free-stream Mach numbers
beta_D = sqrt(1-(M_inf_D.^2));  % Prandtl-Glauert Correction Factor
b_eiw_D = b_0.*beta_D;             % Span of Equivalent Incompressible Wing
Sw_eiw_D = Sw.*beta_D;             % Surface of Equivalent Incompressible Wing

b_kink_eiw_D = zeros(1, length(M_inf_D));

colors = lines(length(M_inf_D)); % At the start, outside the loop
handles = zeros(1, length(M_inf_D)); % Preallocate for legend

for k=1:length(M_inf_D)

% only compress y-coordinates
b = b_eiw_D(k);             % Wingspan [m]
c_root = b_0/3.5;        % Chord @ root [m]
% b_kink = c_root * 1.6 ;   % Span @ kink [m]

b_kink = (1.6/3.5)*b;    % Span @ kink - f(b) not f(c_root)
c_tip = c_root * 0.2;    % Chord @ tip [m]
lambda = 28;             % Swept angle [deg]
lambda_eiw = atan(tand(lambda)/beta_D(k)) * (180/pi); % Equivalent sweep angle in degrees


% Coordinates of characteristic points of the wing
tip = [0,0] ;
b_kink_top = [(b_kink/2) * tand(lambda_eiw) , b_kink/2];
b_kink_bottom = [c_root , b_kink/2];
c_tip_top = [(b/2)*tand(lambda_eiw) , b/2];
c_tip_bottom = [(b/2)*tand(lambda_eiw) + c_tip , b/2];
c_root_profile = [c_root , 0];

% Define the 5 points of the wing outline - 3 sections where geometry changes: Tip, Kink & Root
x_wing_profile = [tip(1), b_kink_top(1), c_tip_top(1), c_tip_bottom(1), b_kink_bottom(1), c_root_profile(1)]; %, tip(1)
y_wing_profile = [tip(2), b_kink_top(2), c_tip_top(2), c_tip_bottom(2), b_kink_bottom(2), c_root_profile(2)]; %, tip(2)
y_wing_profile_left = -y_wing_profile;


%%% Wing geometry and panel generation%%
% Storing panels points in a matrix, the same matrix for all the wing
PanelPoints = zeros((Ny/2)+1, Nx+1, 2); 
PanelPoints_full = zeros(Ny+ 1, Nx+1, 2) ;

% Defining the leading edge coordinates
Leading_edge_root_kink = [tip(1), tip(2); b_kink_top(1), b_kink_top(2)];         % LE Root-Kink
Leading_edge_kink_tip = [b_kink_top(1), b_kink_top(2); c_tip_top(1), c_tip_top(2)]; % LE Kink-Tip

% Defining the trailing edge coordinates
Trailing_edge_root_kink = [c_root_profile(1), c_root_profile(2); b_kink_bottom(1), b_kink_bottom(2)]; % TE Root-Kink
Trailing_edge_kink_tip = [b_kink_bottom(1), b_kink_bottom(2); c_tip_bottom(1), c_tip_bottom(2)];   % TE Kink-Tip

dy = b_eiw_D(k)/Ny;               % Width of each panel
dx = linspace(0, 1, Nx+1); % Length of each panel

for j = 1:(Ny/2) + 1 
    y_j = (j-1) * dy; 

    if y_j > b/2
        y_j = b/2;
    end
    % local_x_le = 0.0;
    % local_x_te = 0.0;

    % Interpolate LE and TE x-coordinates depending if we are before or after the kink
    if y_j <= b_kink_top(2) % Root-kink
        local_x_le = interp1(Leading_edge_root_kink(:,2), Leading_edge_root_kink(:,1), y_j);  % LE x-coordinate root-kink
        local_x_te = interp1(Trailing_edge_root_kink(:,2), Trailing_edge_root_kink(:,1), y_j); % TE x-coordinate
    else % Kink-tip
        local_x_le = interp1(Leading_edge_kink_tip(:,2), Leading_edge_kink_tip(:,1), y_j);  % LE x-coordinate kink-tip
        local_x_te = interp1(Trailing_edge_kink_tip(:,2), Trailing_edge_kink_tip(:,1), y_j); % TE x-coordinate kink-tip
    end

    % Calculates the chord length of each panel
    local_chord = local_x_te - local_x_le;

    % Creates all chordwise panel vertices a this y-coordinate
    for i = 1:(Nx+1)
        PanelPoints(j, i, 1) = local_x_le + dx(i) * local_chord; % x-coordinate of panel i,j
        PanelPoints(j, i, 2) = y_j;         % y-coordinate of panel i,j
    end
end

% Mirroring half of the wing to get the complete wing plot
PanelPoints_full(1:Ny/2+1, :, 2) = -PanelPoints(end:-1:1,:,2) ;
PanelPoints_full(Ny/2+2 : Ny + 1, :, 2) = PanelPoints(2:end,:,2);
PanelPoints_full(1:Ny/2+1, :, 1) = PanelPoints(end:-1:1,:,1);
PanelPoints_full(Ny/2+2 : Ny + 1, :, 1) = PanelPoints(2:end,:,1) ;


% figure; % Plotting the wing surface

col = colors(k, :); % Pick a color for this Mach number

h = plot(y_wing_profile, -x_wing_profile, 'Color', col, 'LineWidth', 1.5);
hold on
plot(y_wing_profile_left, -x_wing_profile, 'Color', col, 'LineWidth', 1.5);
hold on

handles(k) = h; % Store handle for legend

axis equal % Centering axis
hold on

end

hold off
title('Equivalent Incompressible Wings for $M_{\infty}=[0.5-0.95]$', Interpreter='latex', FontWeight='bold', FontSize=14)
legendStrings = arrayfun(@(m) sprintf('M{\\infty} = %.2f', m), M_inf_D, 'UniformOutput', false);
legend(handles, legendStrings, Interpreter= 'latex');
xlabel('y')
ylabel('x')
axis equal
grid on
yticks = get(gca, 'YTick');          % This grabs MATLAB's current y-axis ticks
yticklabels = -yticks;               % Make all tick labels positive
set(gca, 'YTickLabel', yticklabels); % This replaces label -3 with 3, etc.


% PANEL DIVISION FOR EIW
b = sqrt(AR*Sw);     % Wingspan [m]

for k = 1:length(M_inf_D)

beta_D(k) = sqrt(1-(M_inf_D(k)^2));  % Prandtl-Glauert Correction Factor
b_eiw_D (k) =  b_0*beta_D(k);          % Span of Equivalent Incompressible Wing

%hay que repetir todo el proceso del project1 para el equivalent incompressible wing

% only compress y-coordinates
b = b_eiw_D(k);             % Wingspan [m]
c_root = b_0/3.5;        % Chord @ root [m]
% b_kink = c_root * 1.6 ;   % Span @ kink [m]
b_kink = (1.6/3.5)*b; % Span @ kink - f(b) not f(c_root)
c_tip = c_root * 0.2;    % Chord @ tip [m]
lambda = 28;             % Swept angle [deg]
lambda_eiw = atan(tand(lambda)/beta_D(k)) * (180/pi); % Equivalent sweep angle in degrees


% Coordinates of characteristic points of the wing
tip = [0,0] ;
b_kink_top = [(b_kink/2) * tand(lambda_eiw) , b_kink/2];
b_kink_bottom = [c_root , b_kink/2];
c_tip_top = [(b/2)*tand(lambda_eiw) , b/2];
c_tip_bottom = [(b/2)*tand(lambda_eiw) + c_tip , b/2];
c_root_profile = [c_root , 0];

% Define the 6 points of the wing outline - 3 sections where geometry changes: Tip, Kink & Root
x_wing_profile = [tip(1), b_kink_top(1), c_tip_top(1), c_tip_bottom(1), b_kink_bottom(1), c_root_profile(1)]; %, tip(1)
y_wing_profile = [tip(2), b_kink_top(2), c_tip_top(2), c_tip_bottom(2), b_kink_bottom(2), c_root_profile(2)]; %, tip(2)
y_wing_profile_left = -y_wing_profile;


%%% Wing geometry and panel generation %%
% Storing panels points in a matrix, the same matrix for all the wing
PanelPoints = zeros((Ny/2)+1, Nx+1, 2); 
PanelPoints_full = zeros(Ny+ 1, Nx+1, 2) ;

% Defining the leading edge coordinates
Leading_edge_root_kink = [tip(1), tip(2); b_kink_top(1), b_kink_top(2)];         % LE Root-Kink
Leading_edge_kink_tip = [b_kink_top(1), b_kink_top(2); c_tip_top(1), c_tip_top(2)]; % LE Kink-Tip

% Defining the trailing edge coordinates
Trailing_edge_root_kink = [c_root_profile(1), c_root_profile(2); b_kink_bottom(1), b_kink_bottom(2)]; % TE Root-Kink
Trailing_edge_kink_tip = [b_kink_bottom(1), b_kink_bottom(2); c_tip_bottom(1), c_tip_bottom(2)];   % TE Kink-Tip

dy = b_eiw_D(k)/Ny;               % Width of each panel
dx = linspace(0, 1, Nx+1); % Length of each panel

for j = 1:(Ny/2) + 1 
    y_j = (j-1) * dy; 

    if y_j > b/2
        y_j = b/2;
    end
    % local_x_le = 0.0;
    % local_x_te = 0.0;

    % Interpolate LE and TE x-coordinates depending if we are before or after the kink
    if y_j <= b_kink_top(2) % Root-kink
        local_x_le = interp1(Leading_edge_root_kink(:,2), Leading_edge_root_kink(:,1), y_j);  % LE x-coordinate root-kink
        local_x_te = interp1(Trailing_edge_root_kink(:,2), Trailing_edge_root_kink(:,1), y_j); % TE x-coordinate
    else % Kink-tip
        local_x_le = interp1(Leading_edge_kink_tip(:,2), Leading_edge_kink_tip(:,1), y_j);  % LE x-coordinate kink-tip
        local_x_te = interp1(Trailing_edge_kink_tip(:,2), Trailing_edge_kink_tip(:,1), y_j); % TE x-coordinate kink-tip
    end

    % Calculates the chord length of each panel
    local_chord = local_x_te - local_x_le;

    % Creates all chordwise panel vertices a this y-coordinate
    for i = 1:(Nx+1)
        PanelPoints(j, i, 1) = local_x_le + dx(i) * local_chord; % x-coordinate of panel i,j
        PanelPoints(j, i, 2) = y_j;         % y-coordinate of panel i,j
    end
end

% Mirroring half of the wing to get the complete wing plot
PanelPoints_full(1:Ny/2+1, :, 2) = -PanelPoints(end:-1:1,:,2) ;
PanelPoints_full(Ny/2+2 : Ny + 1, :, 2) = PanelPoints(2:end,:,2);
PanelPoints_full(1:Ny/2+1, :, 1) = PanelPoints(end:-1:1,:,1);
PanelPoints_full(Ny/2+2 : Ny + 1, :, 1) = PanelPoints(2:end,:,1) ;


figure; % Plotting the wing surface

col = colors(k, :); % Pick a color for this Mach number

h = plot(y_wing_profile, -x_wing_profile, 'Color', col, 'LineWidth', 1.5);
hold on
plot(y_wing_profile_left, -x_wing_profile, 'Color', col, 'LineWidth', 1.5);
hold on

handles(k) = h; % Store handle for legend

% Plotting panels division inside the wing
for j = 1:Ny+1
    plot(PanelPoints_full(j, :, 2), -PanelPoints_full(j, :, 1), 'b', 'LineWidth', 0.5);
    hold on
end
for i = 1:Nx+1
    plot(PanelPoints_full(:, i, 2), -PanelPoints_full(:, i, 1), 'b', 'LineWidth', 0.5);
    hold on
end
axis equal % Centering axis
hold on

%%% Computing and plotting Control Points
x_c = zeros(Nx, Ny); % Control point x-coordinates
y_c = zeros(Nx, Ny); % Control point y-coordinates

for i = 1:Nx 
    for j = 1:Ny 
        % x,y coordinates for the 4 corners of the (i,j) panel
        P1_coords = squeeze(PanelPoints_full(j,   i,   :)); % Inboard, LE
        P2_coords = squeeze(PanelPoints_full(j,   i+1, :)); % Inboard, TE
        P3_coords = squeeze(PanelPoints_full(j+1, i+1, :)); % Outboard, TE
        P4_coords = squeeze(PanelPoints_full(j+1, i,   :)); % Outboard, LE

        % control points @ exact center of the panel (same as in exercises of sheet 4)
        x_c(i,j) = (P1_coords(1) + P2_coords(1) + P3_coords(1) + P4_coords(1)) / 4;
        y_c(i,j) = (P1_coords(2) + P2_coords(2) + P3_coords(2) + P4_coords(2)) / 4;

        plot(y_c(i,j),-x_c(i,j),'r.', 'MarkerSize', 6) % Adding control points to plot of wing with panels division

    end
end

title(sprintf('Equivalent Incompressible Wing for $M_{\\infty} = %.2f$', M_inf_D(k)), Interpreter='latex', FontWeight='bold',FontSize=14)
legendStrings = arrayfun(@(m) sprintf('M{\\infty} = %.2f', m), M_inf_D, 'UniformOutput', false);
legend(handles, legendStrings, Interpreter= 'latex');
xlabel('y (spanwise direction)', Interpreter= 'latex')
ylabel('x (chordwise direction)', Interpreter= 'latex')
axis on
grid off
limits = [min(y_wing_profile_left)-1,max(y_wing_profile)+1];
xlim (limits)
yticks = get(gca, 'YTick');      % This grabs MATLAB's current y-axis ticks
yticklabels = -yticks;       % Make all tick labels positive
set(gca, 'YTickLabel', yticklabels); % This replaces label -3 with 3, etc.
hold off

% end

end

%% VLM FOR EIW TO GET dCL/dα - from Project 1 - OK revisado con Ignacio - now for M_inf_D=[0.5-0.95]
alpha_deg = [0, 5, 10, 15, 20] ; % Angle of attack [degrees]
alpha_rad = deg2rad(alpha_deg);  % Angle of attack [radians]

Na = length(alpha_rad);
U_inf_D = zeros(1,length(M_inf_D));

% Preallocating vectors for Lift
CL_i_D = zeros(Na);
CL_c_D = zeros(Na);
L_D = zeros(Na);
CL_full_i_D = zeros(Na,Na);
CL_full_c_D = zeros(Na,Na);
dCL_dalpha_c_rad_D = zeros(1,Na);
dCL_dalpha_c_deg_D = zeros(1,Na);

% Preallocating vectors for plot
colors = lines(length(M_inf_D));     % Vector for different colors for each Mach
handles = zeros(1, length(M_inf_D)); % Vector for legend

% Preallocating vectors for Drag
D0 = zeros(1, length(M_inf_D));         % D0 ≡ Parasitic Drag
% D0_PG = zeros(1, length(M_inf_D));    % D0_PG ≡ Parasitic Drag w/ PG contraction
Di_i_y = zeros(Ny, length(M_inf_D));    % Di_i_y ≡ Induced drag in each panel
Di_i = zeros(1, length(M_inf_D));       % Di_i ≡ Induced Drag of Equivalent Incompressible Wing
Di_c = zeros(1, length(M_inf_D));       % Di_i ≡ Induced Drag of Equivalent Incompressible Wing
D_total = zeros(1, length(M_inf_D));    % D_total ≡ Total drag = Parasitic Drag + Induced drag
q_inf = zeros(1, length(M_inf_D));      % q_inf ≡ Dynamic Pressure

% Preallocating vectors for required AoA to maintain steady level flight
alpha_req = zeros(1,length(M_inf_D));
CL_req = zeros(length(M_inf_D),1);

% Coordinates of Trefftz's Plane - only y used
y_trefftz = zeros(Ny+1,length(M_inf_D));
y_trefftz_cp = zeros(Ny,length(M_inf_D));
dy_i = zeros(1,length(M_inf_D));

% Preallocating vectors for Gamma & w
Gamma_D_full = zeros (Nx*Ny,length(M_inf_D));
Gamma_trefftz = zeros(Ny,length(M_inf_D)); % (span, Mach)
Gamma_trefftz_delta = zeros(Ny+1,length(M_inf_D)); % (span, Mach)
w_trefftz = zeros(Ny,length(M_inf_D));

for i = 1:length(M_inf_D)

        dy_i (i) = b_eiw_D(i)/Ny;

        % LIFT CALCULATION 
        [CL_i_D, ww_mat_D, Gamma_D, L_D, U_inf_d] = vlm (M_inf_D(i), PanelPoints_full, alpha_rad, Nx, Ny, x_c, y_c, b_eiw_D(i), rho, a, Sw_eiw_D(i));

        U_inf_D (i) = U_inf_d;
     
        CL_full_i_D(i,:) = CL_i_D;    % Store all values of CL - Each row is for each M_inf_D   -->  Equal for all M_inf_D
        CL_c_D = CL_i_D./beta_D(i);     % Undo PG contraction for CL
        CL_full_c_D(i,:) = CL_c_D;    % Store all values of CL - Each row is for each M_inf_D

        col2 = colors(i, :);          % Pick a color for plotting CL for this Mach number - same colors as EIWs
        h2 = plot(alpha_deg, CL_c_D, 'Color', col2, 'LineWidth', 1.5);
        hold on
        handles(i) = h2;          % Store handles for legend

        dCL_dalpha_c_rad_D(i) = (CL_full_c_D(i, end) - CL_full_c_D(i, 1)) / (alpha_rad(end) - alpha_rad(1)); 
        dCL_dalpha_c_deg_D(i) = (CL_full_c_D(i, end) - CL_full_c_D(i, 1)) / (alpha_deg(end) - alpha_deg(1));

        % PARASITIC DRAG
        q_inf(i) = 1/2 * rho * U_inf_D(i)^2;   % Dynamic pressure
        D0(i) = (CD0 * q_inf(i) * Sw);  % Parasitic drag
        % D0_PG(i) = (CD0 / beta_D(i) * q_inf(i) * Sw_eiw_D(i));  % Parasitic drag w/ PG contraction -> Sw_eiw_D & CD(c) = CD(i)/beta_D
        % since beta is multiplying CD0 and dividing Sw - D0_PG≡D0
       
        % ANGLE OF ATTACK REQUIRED TO MAINTAIN STEADY LEVEL FLIGHT

        % CL required for steady levelled flight W=L=CL*q_inf*Sw
        CL_req (i) = W / (q_inf(i) * Sw);
        alpha_req(i) = interp1(CL_full_c_D(i,:), alpha_deg, CL_req(i), 'linear', 'extrap'); % interpolating to find alpha such that CL(alpha) = CL_req(i)

        % SUGGESTION BY PERPLEXITY - COMPUTE GAMMA ONLY FOR THE REQUIRED AoA
        alpha_req_rad = deg2rad(alpha_req(i));                      % NEW

        % RHS for induced w at required alpha                       % NEW
        b_w_req = -U_inf_D(i) * sin(alpha_req_rad) * ones(Nx*Ny,1); % NEW

        % Solve for Gamma at required alpha                         % NEW
        Gamma_req = ww_mat_D \ b_w_req;                             % NEW

        % Store and collapse chordwise to Trefftz distribution      % NEW
        Gamma_D_full(:, i) = Gamma_req;                             % NEW
        Gamma_mat_req = reshape(Gamma_req, Ny, Nx);                 % NEW
        Gamma_trefftz(:, i) = sum(Gamma_mat_req, 2);                % NEW
        
        % % Computing Gamma for Trefftz's Plane
        % Gamma_D_full (:,i) = Gamma_D(:);        % Trailing-edge circulation for each spanwise location - needed for Trefftz
        % 
        % % Let Gamma_D_full be (Ny*Nx, Nmach)
        % 
        % Gamma_vec = Gamma_D_full(:, i); % (Ny*Nx) x 1 vector for Mach m
        % Gamma_mat = reshape(Gamma_vec, Ny, Nx); % Ny x Nx, spanwise major
        % Gamma_trefftz(:, i) = sum(Gamma_mat, 2); % sum along columns (chordwise stations)

        % INDUCED DRAG CALCULATION    
        % computing dGamma for the integral of w
        for j = 1:Ny+1
            if j==1
                Gamma_trefftz_delta(j,i) = Gamma_trefftz(j,i);
            elseif j == (Ny+1)
                Gamma_trefftz_delta(j,i) = - Gamma_trefftz(j-1,i);
            else
                Gamma_trefftz_delta(j,i) = Gamma_trefftz(j,i) - Gamma_trefftz(j-1,i);
            end
        end

        % Computing position of control points in the middle of each Trefftz's Plane panel
        y_trefftz (:,i) = linspace(-b_eiw_D(i)/2, b_eiw_D(i)/2, Ny+1);  % Contracting span depending on Mach
        y_trefftz_cp = 0.5 * (y_trefftz(1:Ny, :) + y_trefftz(2:Ny+1, :));
        % for j = 1:Ny - same procedure
        %     y_trefftz_cp (j,i) = 1/2 * (y_trefftz(j+1,i)-y_trefftz(j,i)) + y_trefftz(j);
        % end

        for k = 1:Ny
            w_trefftz(k,i) = 0; % avoid accumulation errors
            for j = 1:Ny+1
                   w_trefftz(k,i) = -1/(2*pi) * Gamma_trefftz_delta(j,i)/((y_trefftz_cp(k,i))-y_trefftz(j,i)) + w_trefftz(k,i);
            end
        end

        % Solving integral to get Induced Drag
        for j = 1:Ny
            Di_i(i) = -1/2 * rho * Gamma_trefftz(j,i) * w_trefftz(j,i) * dy_i(i) + Di_i(i);
        end

        Di_c(i) = Di_i(i)/beta_D(i);      % Induced drag for actual compressible regime - undoing PG contraction
        D_total(i) = D0(i) + Di_c(i);     % Total drag
end

% Plot of CL vs alpha
title('$C_L$ vs $\alpha$ for $M_{\infty}=[0.5-0.95]$', Interpreter='latex', FontWeight='bold', FontSize=14)
legendStrings = arrayfun(@(m) sprintf('M{\\infty} = %.2f', m), M_inf_D, 'UniformOutput', false);
legend(handles, legendStrings, Interpreter= 'latex', Location='best');
xlabel('$\alpha$', Interpreter='latex', FontSize=14)
ylabel('$C_L$', Interpreter='latex', FontSize=14)
grid on
hold off

% Plotting ALL contributions of DRAG
% Plot of D0, Di_c & D_total vs Uinf
figure
plot(U_inf_D, D0, 'b','LineWidth', 2.5, 'DisplayName', 'Parasitic drag - $D_0$');         % Parasitic Drag vs. Free-stream velocity
hold on
plot(U_inf_D, Di_c, 'r','LineWidth', 2.5, 'DisplayName', 'Induced Drag - $D_i$');         % Induced Drag vs. Free-stream velocity
hold on
plot(U_inf_D, D_total, 'k','LineWidth', 2.5, 'DisplayName', 'Total Drag - $D_0 + D_i$');  % Total Drag vs. Free-stream velocity
legend('Location', 'best', 'FontSize', 12, 'FontWeight', 'normal', 'Interpreter','latex');
ylabel('$D [N]$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$U_\infty [m/s]$', 'Interpreter', 'latex', 'FontSize', 12)
title('Drag forces vs. Free-stream velocity','FontSize', 12, 'Interpreter', 'latex');
hold off
grid on;

%% Computation of optimal flight speed (U_opt) - U_inf that results in minimum total

% [minDrag, index] = min(D_total);    % Finding minimum value of Total Drag & its position
% U_opt_ms = U_inf_D(index);                % Optimal flight speed [m/s]
% % U_opt_ms2 = interp1(CL_full_c_D(i,:), alpha_deg, CL_req(i), 'linear', 'extrap'); % interpolating to find min Drag
% U_opt_kmh = U_opt_ms * 3600/1000;   % Optimal flight speed [km/h] --> km/h = m/s * s/h * km/m

% Smooth interpolation of D_total(U_inf_D)
D_interp = @(U) interp1(U_inf_D, D_total, U, 'spline');
% Search interval: between your min and max speeds
Umin = min(U_inf_D);
Umax = max(U_inf_D);
% Find continuous minimum
U_opt_ms = fminbnd(D_interp, Umin, Umax);
D_opt = D_interp(U_opt_ms);
U_opt_kmh = U_opt_ms * 3.6;

% REQUIRED AoA FOR OPTIMAL FLIGHT SPEED
alpha_opt_req_deg = interp1(U_inf_D, alpha_req, U_opt_ms, 'spline');
alpha_opt_req_rad = deg2rad(alpha_opt_req_deg);

% Plotting Optimal Flying Speed in Drag components vs. Free-stream velocity Figure
figure;
plot(U_inf_D, D0, 'b', 'LineWidth', 2.5, 'DisplayName', 'Parasitic Drag - $D_0$');            % Parasitic Drag vs. Free-stream velocity    
hold on;
plot(U_inf_D, Di_c, 'r', 'LineWidth', 2.5, 'DisplayName', 'Induced Drag - $D_i$');              % Induced Drag vs. Free-stream velocity
hold on
plot(U_inf_D, D_total, 'k', 'LineWidth', 2.5, 'DisplayName', 'Total Drag - $D_0 + D_i$');     % Total Drag vs. Free-stream velocity 
hold on
xline(U_opt_ms, '--', 'LineWidth', 1.5, 'DisplayName', 'Optimal Flight Speed');        % Vertical line @ Optimal Flight Speed
hold on
yline(D_opt, '--g', 'LineWidth', 1.5, 'DisplayName', 'Minimum Total Drag');           % Horizontal line @ Minimum Total Drag
hold on
plot(U_opt_ms, D_opt, 'ko', 'MarkerSize', 11, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');        % Intersection of min.Drag & Opt.Flight Speed
text(U_opt_ms + 5, D_opt - 4000, sprintf('$U_{opt} = %.2f \\; \\mathrm{m/s}$\n$D_{min} = %.2f \\; \\mathrm{N}$', U_opt_ms, D_opt), ...
    'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'latex', 'HorizontalAlignment', 'left');          % Values of Min.Drag & Opt.Flying Speed
legend('Location', 'best', 'FontSize', 12, 'FontWeight', 'normal', 'Interpreter','latex');
ylabel('$D [N]$', 'Interpreter', 'latex');
xlabel('$U_{\infty} [m/s]$', 'Interpreter', 'latex');
title('Optimal Flight Speed for $h=11km$', 'Interpreter', 'latex', FontSize=14);
hold off;
grid on;
%% Aerodynamic Efficiency

L = W;                 % Lift @ Steady Levelled Flight
% E_0 = L./D0;         % Aerodynamic efficiency for Parasitic Drag = Lift/ Parasitic Drag
% E_i = L./Di_c;         % Aerodynamic efficiency for Induced Drag = Lift/ Induced Drag
% E_total = L./D_total;  % Aerodynamic efficiency for Total Drag = Lift/ Total Drag

D0_opt   = interp1(U_inf_D, D0,   U_opt_ms, 'spline');
Di_opt   = interp1(U_inf_D, Di_c, U_opt_ms, 'spline');
Dtot_opt = interp1(U_inf_D, D_total, U_opt_ms, 'spline');   % or just minDrag

E_0     = L / D0_opt;
E_i     = L / Di_opt;
E_total = L / Dtot_opt;

% figure
% plot(U_inf_D,E_0, 'b','LineWidth', 2.5, 'DisplayName', '$E(D_0)$');                 % Aerodynamic Efficiency for Parasitic Drag vs. Free-stream velocity
% hold on
% plot(U_inf_D, E_i, 'y','LineWidth', 2.5, 'DisplayName', '$E(D_i)$');                % Aerodynamic Efficiency for Induced Drag vs. Free-stream velocity
% hold on
% plot(U_inf_D, E_total, 'r','LineWidth', 2.5, 'DisplayName', '$E(D_0 + D_i)$')       % Aerodynamic Efficiency for Total Drag vs. Free-stream velocity
% legend('Location', 'best', 'FontSize', 12, 'FontWeight', 'normal', 'Interpreter','latex');
% ylabel('$E = C_L/C_D \ [-]$', 'Interpreter', 'latex');
% xlabel('$U_\infty [m/s]$', 'Interpreter', 'latex')
% title('Aerodynamic efficiencies vs. Free-stream velocity', 'Interpreter', 'latex', 'FontWeight', 'bold','FontSize', 13);
% hold off
% grid on;

%% Required AoA for steady levelled flight

% alpha_req = zeros(1,length(M_inf_D));
% CL_req = zeros(length(M_inf_D),1);
% CL_diff = zeros(size(CL_full_c_D));
% min_CL_diff = zeros(length(M_inf_D));
% 
% for i = 1:length(M_inf_D)
%     % CL required for steady levelled flight W=L=CL*q_inf*Sw
%     CL_req (i) = W / (q_inf(i) * Sw);
%     alpha_req(i) = interp1(CL_full_c_D(i,:), alpha_deg, CL_req(i), 'linear', 'extrap'); % interpolating to find alpha such that CL(alpha) = CL_req(i)
% end

figure;
plot(U_inf_D, alpha_req,'LineWidth',2, 'Marker','o', 'Color',"r")
xlabel('$U_\infty$ $[m/s]$', 'Interpreter', 'latex');
ylabel('Angle of Attack, $\alpha$ [deg]', 'Interpreter', 'latex');
title('Required Angle of Attack to maintain steady level flight at $h=11km$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
grid on

%% Functions

function [uu,vv,ww] = vring(Gamma,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,xx,yy,zz) 

%function [uu,vv,ww] = vring(Gamma,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,xx,yy,zz)
% Computes the velocities produced in (xx,yy,zz) by a QUAD vortex ring 
% of circulation Gamma with vertices (x,y,z)_{1,2,3,4}. 
%
% Based on the formulas in Katz's book, page 251-256

x0 = [x1 x2 x3 x4 x1]; 
y0 = [y1 y2 y3 y4 y1]; 
z0 = [z1 z2 z3 z4 z1]; 

tol = 1d-10; 

uu = zeros(size(xx)); 
vv = zeros(size(xx)); 
ww = zeros(size(xx)); 

for i=1:4

    curlx = (yy-y0(i)).*(zz-z0(i+1)) - (zz - z0(i)).*(yy-y0(i+1)); 
    curly =-(xx-x0(i)).*(zz-z0(i+1)) + (zz - z0(i)).*(xx-x0(i+1)); 
    curlz = (xx-x0(i)).*(yy-y0(i+1)) - (yy - y0(i)).*(xx-x0(i+1)); 
    rr1   = sqrt((xx-x0(i)).^2 + (yy-y0(i)).^2 + (zz-z0(i)).^2); 
    rr2   = sqrt((xx-x0(i+1)).^2 + (yy-y0(i+1)).^2 + (zz-z0(i+1)).^2); 
    dot1  = (x0(i+1)-x0(i)).*(xx-x0(i  )) + (y0(i+1)-y0(i)).*(yy-y0(i  )) + (z0(i+1)-z0(i)).*(zz-z0(i  )); 
    dot2  = (x0(i+1)-x0(i)).*(xx-x0(i+1)) + (y0(i+1)-y0(i)).*(yy-y0(i+1)) + (z0(i+1)-z0(i)).*(zz-z0(i+1)); 

    curl2 = (curlx.^2 + curly.^2 + curlz.^2);

    KK    = Gamma./(4*pi.*(curlx.^2 + curly.^2 + curlz.^2)).*(dot1./rr1 - dot2./rr2); 
    KK(find(rr1<tol | rr2<tol | curl2<tol)) = 0; 

    uu    = uu + KK.*curlx; 
    vv    = vv + KK.*curly; 
    ww    = ww + KK.*curlz; 

end


return
end

function [CL_i, ww_mat, Gamma, L, V_inf] = vlm (M_inf, PanelPoints_full, alpha_rad, Nx, Ny, x_c, y_c, b_eiw, rho, a, Sw_eiw)
    
    if M_inf == 0 % just to avoid having V=a*(M=0)
        V_inf = 1;
    else
        V_inf = M_inf*a;
    end

    ww_mat = zeros(Nx, Ny);
    idx_cp=0;

    for m = 1:Nx
        for n = 1:Ny


        idx_cp = idx_cp+1; %(n - 1) * Nx + m;  % Control point index

        % Control point coordinates
        xx_cp = x_c(m, n);
        yy_cp = y_c(m, n);
        zz_cp = 0;

        idx_vr =0;
        for i = 1:Nx
            for j = 1:Ny

                idx_vr = idx_vr+1;  % Vortex ring index

                Gamma = 1;

                % Panel corner coordinates
                x1 = PanelPoints_full(j, i, 1);
                y1 = PanelPoints_full(j, i, 2);
                z1 = 0;

                x2 = 10000;
                %PanelPoints_full(j, i+1, 1);
                y2 = PanelPoints_full(j, i+1, 2);
                z2 = 0;

                x3 = 10000;
                %PanelPoints_full(j+1, i+1, 1);
                y3 = PanelPoints_full(j+1, i+1, 2);
                z3 = 0;

                x4 = PanelPoints_full(j+1, i, 1);
                y4 = PanelPoints_full(j+1, i, 2);
                z4 = 0;



                xx = x_c(i,j);
                yy = y_c(i,j);
                zz = 0;


                [uu, vv, ww] = vring(Gamma, x1, x4, x3, x2, y1, y4, y3, y2, z1, z2, z3, z4, xx_cp, yy_cp, zz_cp);



                ww_mat(idx_cp, idx_vr) =  ww;

                % L(i,j) = rho*induced_velocity(i,j,3);

            end
        end
        end
    end


    % Computing w & Gamma

    for k = 1:length(alpha_rad)
    alpha = alpha_rad(k);

    b_w = -V_inf * sin(alpha) * ones(Nx * Ny, 1);   % ~alpha??

    % Solve the linear system for Gamma
    Gamma = ww_mat \ b_w;


    dy_i = (b_eiw) / Ny;  % spanwise panel width
    L(k) = rho * V_inf * sum(Gamma) * dy_i;
    % L(k) = sum(Lp);
    CL_i(k) = L(k) / (0.5 * rho * V_inf^2 * Sw_eiw);
    end
end