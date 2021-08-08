%% Data form

% INPUT DATA
% x = x-coordinate of trajectory
% y = y-coordinate of trajectory
% name = 'name_of_track'
% m = mass
% ftmax = max traction (N)
% fbmax = max braking (-ve)(N)
% fnmax = max cornering (N)

% OUTPUT DATA
% velProf = velocity profile of the given track

function [velProf] = func_velProf(x,y,name,m,ftmax,fbmax,fnmax)
%% Initialization

n = numel(x);
accel = zeros(size(x)); % acceleration profile
decel = zeros(size(x)); % deceleration profile
% ftmax = ftmax*m; % traction max
% fbmax = -fbmax*m; % braking max
% fnmax = fnmax*m; % cornering max
drag = 0.0021*m; % drag

%% Segment length

len = zeros(size(x));

for i = 2:n
    len(i) = len(i-1)+sqrt((x(i)-x(i-1))^2+(y(i)-y(i-1))^2);
end

%% Curvature

[~,R,~] = curvature([x y]);
K = 1./R;
K(1) = 0;
K(end) = 0;
[~,locs] = findpeaks(K);

%% Velocity Profile Calculation

% start section
accel(1) = 0;
for i = 2:locs(1)
     accel(i) = sqrt(2*sqrt((ftmax-drag*accel(i-1)^2)^2 ...
                            -(m*K(i-1)*((ftmax-drag*accel(i-1)^2)/fnmax)^2*accel(i-1)^2)^2)*(len(i)-len(i-1))/m ...
                            + accel(i-1)^2);
end
    
for i = 1:numel(locs)
    % check if accel vel is greater than crit vel
    vcrit = sqrt(fnmax/(m*K(locs(i))));
    if accel(locs(i)) < vcrit
        curVel = accel(locs(i));
    else
        curVel = vcrit;
        accel(locs(i)) = vcrit;
    end

    % next step using cur vel
    accel(locs(i)+1) = sqrt(2*sqrt((ftmax-drag*curVel^2)^2 ...
                                    -(m*K(locs(i)-1)*((ftmax-drag*curVel^2)/fnmax)^2*curVel^2)^2)*(len(locs(i))-len(locs(i)-1))/m ...
                                    + curVel^2);
    
    if i == numel(locs)
        % end section
        for j = locs(i)+2:n
            accel(j) = sqrt(2*sqrt((ftmax-drag*accel(i-1)^2)^2 ...
                                    -(m*K(j-1)*((ftmax-drag*accel(i-1)^2)/fnmax)^2*accel(j-1)^2)^2)*(len(j)-len(j-1))/m ...
                                    + accel(j-1)^2);
        end
    else
        % usual step
        for j = locs(i)+2:locs(i+1)
            accel(j) = sqrt(2*sqrt((ftmax-drag*accel(j-1)^2)^2 ...
                                    -(m*K(j-1)*((ftmax-drag*accel(j-1)^2)/fnmax)^2*accel(j-1)^2)^2)*(len(j)-len(j-1))/m ...
                                    + accel(j-1)^2);
        end
    end
end

% end point
decel(end) = 200; % max final velocity 
for j = numel(decel)-1:-1:locs(end)
    decel(j) = sqrt(2*sqrt((fbmax-drag*decel(j+1)^2)^2 ...
                            -(m*K(j+1)*((fbmax-drag*decel(j+1)^2)/fnmax)^2*decel(j+1)^2)^2)*(len(j+1)-len(j))/m ...
                            + decel(j+1)^2);
end

for i = numel(locs):-1:1
    % check if decel vel is greater than crit vel
    vcrit = sqrt(fnmax/(m*K(locs(i))));
    if decel(locs(i)) < vcrit
        curVel = decel(locs(i));
    else
        curVel = vcrit;
        decel(locs(i)) = vcrit;
    end
    % next step using cur vel
    decel(locs(i)-1) = sqrt(2*sqrt((fbmax-drag*curVel^2)^2 ...
                                    -(m*K(locs(i))*((fbmax-drag*curVel^2)/fnmax)^2*curVel^2)^2)*(len(locs(i))-len(locs(i)-1))/m ...
                                    + curVel^2);

    if i == 1
        % start section
        for j = locs(i)-2:-1:1
            decel(j) = sqrt(2*sqrt((fbmax-drag*decel(j+1)^2)^2 ...
                                    -(m*K(j+1)*((fbmax-drag*decel(j+1)^2)/fnmax)^2*decel(j+1)^2)^2)*(len(j+1)-len(j))/m ...
                                    + decel(j+1)^2);
        end
    else
        % usual step
        for j = locs(i)-2:-1:locs(i-1)
            decel(j) = sqrt(2*sqrt((fbmax-drag*decel(j+1)^2)^2 ...
                                    -(m*K(j+1)*((fbmax-drag*decel(j+1)^2)/fnmax)^2*decel(j+1)^2)^2)*(len(j+1)-len(j))/m ...
                                    + decel(j+1)^2);
        end
    end
end

velProf = min(accel,decel);

%% Lap Time calcluation

time = zeros(size(x));
for i = 2:numel(x)
    acc = (velProf(i)^2-velProf(i-1)^2)/(2*(len(i)-len(i-1)));
    time(i) = time(i-1) + (velProf(i)-velProf(i-1))/acc;
end

%% Plot Velocity Profile

% figure
% plot(accel)
% hold on
% plot(decel)

figure
hold on
plot(len,velProf*3.6,'LineWidth',2)
grid on
xlabel('s(m)','fontweight','bold','fontsize',14)
ylabel('kmph','fontweight','bold','fontsize',14)
title(sprintf('%s - Velocity Profile\nLap Time = %.2fs',name,time(end)),'fontsize',16)
