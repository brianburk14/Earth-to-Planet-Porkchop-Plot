clear; clc; clf; close all;

%Read the data files and generate the state vectors as well as the dates
%for Earth and Mars at each time
[state_array_EARTH, julian_date_line_EARTH] = stateReader('2005_Earth_Dep.txt');
[state_array_MARS, julian_date_line_MARS] = stateReader('2005_Mars_Arr.txt');

C3_levels = [15.5, 16, 16.5, 17, 18, 20, 22.5, 25, 30];
%C3_levels = 15:4:100;

%Assign a step size in days to get an accurate porkchop plot as well as
%have the program not take too long to run
step_time = 4;
index_EARTH = 0;
for k = 1:julian_date_line_EARTH
    if(k == step_time*index_EARTH + 1)
        state_used_EARTH(index_EARTH+1, :) = state_array_EARTH(k, :);
        index_EARTH = index_EARTH + 1;
    end
end

index_MARS = 0;
for k = 1:julian_date_line_MARS
    if(k == step_time*index_MARS + 1)
        state_used_MARS(index_MARS+1, :) = state_array_MARS(k, :);
        index_MARS = index_MARS + 1;
    end
end

%Execute the Lambert fit, the radial distance from Earth and Mars, as well
%as the transfer time for each case using nested for loop
for i = 1:length(state_used_EARTH)
    for j = 1:length(state_used_MARS)
        transfer_time = (state_used_MARS(j,1)-state_used_EARTH(i,1))*86400;
        if(transfer_time < 0)
            continue;
        end
        [v1, v2] = LambertSolver([state_used_EARTH(i,2); state_used_EARTH(i,3);...
            state_used_EARTH(i,4)], [state_used_MARS(j,2); state_used_MARS(j,3);...
            state_used_MARS(j,4)], transfer_time);
        delta_vD(i,j) = abs(norm(v1 - [state_used_EARTH(i,5); state_used_EARTH(i,6);...
            state_used_EARTH(i,7)]));
        delta_vA(i,j) = abs(norm(v2) - norm([state_used_MARS(j,5); state_used_MARS(j,6);...
            state_used_MARS(j,7)]));
        V1(i, j) = norm(v1);
        C3(i,j) = abs(delta_vD(i,j))^2;
        departure_date(i) = state_used_EARTH(i,1);
        arrival_date(j) = state_used_MARS(j,1);
    end
end
clc;

%Create a mesh grid for the departure and arrival dates
[Y, X] = meshgrid(arrival_date, departure_date);
transfer_matrix = Y - X;

%Convert the julian dates to actual dates
departure_date_value = datestr(datetime(departure_date', "ConvertFrom",...
    "juliandate"), "mm/dd/yyyy");
arrival_date_value = datestr(datetime(arrival_date', "ConvertFrom",...
    "juliandate"), "mm/dd/yyyy");

%Create a mesh grid now using the actual date format
[Y, X] = meshgrid(datenum(arrival_date_value), datenum(departure_date_value));
cmap = jet(length(C3_levels));

[n, m] = size(transfer_matrix);
for i = 1:n
    for j = 1:m
        if(transfer_matrix(i,j) < 0)
            transfer_matrix(i,j) = NaN;
        end
    end
end

%Plot the porkchop plots along with the time of flights
figure(1)
[C, h] = contour(X, Y, C3, C3_levels);
grid on;
hold on;
contour(X, Y, transfer_matrix);
colormap(cmap);

colorbar;
xlabel("Departure Date");
ylabel("Arrival Date");
clabel(C,h);
title("C3 Departure: Earth-Mars");
datetick('x', 'mm/dd/yyyy', 'keeplimits');
datetick('y', 'mm/dd/yyyy', 'keeplimits');
% hold on;
% plot(datenum('08/25/2022'), datenum('07/11/2023'), '-s','MarkerSize',10,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .6]);

function [state_array, julian_date_line] = stateReader(filename)
    start_line = 'SOE';
    end_line = 'EOE';

    N = 0;
    f = fopen(filename,'rt');
    
    julian_date_line = 0;
    orbit_position_line = 0;
    velocity_line = 0;

    while true
        dataline = fgetl(f);
        N = N + 1;

        if (ischar(dataline) == 1)
            if(ischar(dataline) && isempty(regexp(dataline,start_line, 'once')) && isempty(regexp(dataline,end_line, 'once')))
                if(N == 4*julian_date_line + 2)
                    julian_date_line = julian_date_line + 1;
                    julian_cell = split(dataline);
                    julian_date_value = str2double(cell2mat(julian_cell(1)));
                    julian_date_array(julian_date_line) = julian_date_value;
                elseif(N == 4*orbit_position_line + 3)
                    orbit_position_line = orbit_position_line + 1;
                    orbit_position_cell = split(dataline, [" ", "="]);
                    orbit_position_cell = orbit_position_cell(~cellfun('isempty',orbit_position_cell));
                    x_cell = str2double(cell2mat(orbit_position_cell(2)));
                    x_array(orbit_position_line) = x_cell;
    
                    y_cell = str2double(cell2mat(orbit_position_cell(4)));
                    y_array(orbit_position_line) = y_cell;
    
                    z_cell = str2double(cell2mat(orbit_position_cell(6)));
                    z_array(orbit_position_line) = z_cell;
                elseif(N == 4*velocity_line + 4)
                    velocity_line = velocity_line + 1;
                    velocity_cell = split(dataline, [" ", "="]);
                    velocity_cell = velocity_cell(~cellfun('isempty',velocity_cell));

                    vx_cell = str2double(cell2mat(velocity_cell(2)));
                    vx_array(velocity_line) = vx_cell;

                    vy_cell = str2double(cell2mat(velocity_cell(4)));
                    vy_array(velocity_line) = vy_cell;

                    vz_cell = str2double(cell2mat(velocity_cell(6)));
                    vz_array(velocity_line) = vz_cell;
                end
            end

        end

        if (ischar(dataline) == 0)
            break; 
        end 
    end
    fclose(f);
    state_array = [julian_date_array', x_array', y_array', z_array', vx_array', vy_array', vz_array'];
end

function [v1, v2] = LambertSolver(r1, r2, delta_t)
    mu = 1.3271E11;
    r1_cross_r2_z = dot(cross(r1, r2), [0; 0; 1]);

    %Prograde Orbit
    if(r1_cross_r2_z >= 0)
        delta_theta = acos(dot(r1, r2)/(norm(r1)*norm(r2)));
    elseif(r1_cross_r2_z < 0)
        delta_theta = 2*pi - acos(dot(r1, r2)/(norm(r1)*norm(r2)));
    end
    
    A = (sin(delta_theta))*sqrt((norm(r1)*norm(r2))/(1 - cos(delta_theta)));

    tol = 1.e-8;
    nmax = 500;
    n = 0;
    ratio = 1;

    z = -100;
    while F_fun(z, delta_t, r1, r2, A, mu) < 0
        z = z + 0.1;
    end

    while (abs(ratio) > tol) && (n <= nmax)
        n = n + 1;  
        ratio = F_fun(z, delta_t, r1, r2, A, mu)/dF_fun(z, r1, r2, A);
        z = z - ratio;
    end

    f = 1 - y_fun(z, r1, r2, A)/norm(r1);
    g = A*sqrt(y_fun(z, r1, r2, A)/mu);
    g_dot = 1 - y_fun(z, r1, r2, A)/norm(r2);

    v1 = 1/g*(r2 - f*r1);
    v2 = 1/g*(g_dot*r2 - r1);
end

function sol = S_fun(z)
    if(z == 0)
        sol = 1/6;
    elseif(z > 0)
        sol = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
    elseif(z < 0)
        sol = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    end
end

function sol = C_fun(z)
    if(z == 0)
        sol = 1/2;
    elseif(z > 0)
        sol = (1 - cos(sqrt(z)))/z;
    elseif(z < 0)
        sol = (cosh(sqrt(-z)) - 1)/(-z);
    end
end

function sol = y_fun(z, r1, r2, A)
    sol = norm(r1) + norm(r2) + A*(z*S_fun(z) - 1)/sqrt(C_fun(z));
end

function sol = F_fun(z, delta_t, r1, r2, A, mu)
    sol = (y_fun(z, r1, r2, A)/C_fun(z))^(3/2)*S_fun(z) + A*sqrt(y_fun(z, r1, r2, A)) - sqrt(mu)*delta_t;
end

function sol = dF_fun(z, r1, r2, A)
    if(z ~= 0)
        sol = (y_fun(z, r1, r2, A)/C_fun(z))^(3/2)*((C_fun(z) - 3*S_fun(z)/(2*C_fun(z)))/(2*z) +...
            3*S_fun(z)^2/(4*C_fun(z))) + A/8*(3*S_fun(z)/C_fun(z)*sqrt(y_fun(z, r1, r2, A)) + A*sqrt(C_fun(z)/y_fun(z, r1, r2, A)));
    else
        sol = sqrt(2)/40*y_fun(z, r1, r2, A)^(3/2) + A/8*(sqrt(y_fun(z, r1, r2, A)) + A*sqrt(1/(2*y_fun(z, r1, r2, A))));
    end
end