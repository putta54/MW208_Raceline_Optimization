% Install following Add-On for curvaure function to work
% https://in.mathworks.com/matlabcentral/fileexchange/69452-curvature-of-a-1d-curve-in-a-2d-or-3d-space 
%% Data form

% INPUT DATA
% name = 'name_of_track'
% m = mass
% ftmax = max traction (N)
% fbmax = max braking (N)
% fnmax = max cornering (N)
% trackData - [x-ref y-ref xin yin xout yout]
% Values I had taken
% m = 740, ftmax = 16*740, fbmax = -18*740, fnmax = 30*740

% OUTPUT DATA
% velProf = velocity profile of the given track

function [velProf] = func_velProf(traj,name,m,ftmax,fbmax,fnmax,trackData)
%% Initialization

x = traj(:,1);
y= traj(:,2);
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

%% Plot Velocity Profile (v vs s)

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

%% Plot velocity profile onto trajectory

figure
plot(traj(:,1),traj(:,2),'color','w','linew',1)
hold on

% plot starting line
plot([trackData(1,3) trackData(1,5)], [trackData(1,4) trackData(1,6)],'color','b','linew',2)
% plot([xin(2) xout(2)], [yin(2) yout(2)],'color','k','linew',2)

% plot reference line
plot(trackData(:,1),trackData(:,2),'--')
hold on

% plot inner track
plot(trackData(:,3),trackData(:,4),'color','k','Linewidth',0.2)

%plot outer track
plot(trackData(:,5),trackData(:,6),'color','k','Linewidth',0.2)

scatter(traj(:,1),traj(:,2),5,velProf,'filled')
colormap('parula');
colBar = colorbar('eastoutside');
colBar.Label.String = "Velocity [m/s]";
xlabel('X [m]');
ylabel('Y [m]');
% title(['Best Lap Velocity-Curvature Optimization on ',testName],['Iteration: ',num2str(videoSet{end,2}),' | Curvature: k=',num2str(videoSet{end,3}), '[1/m] | Minimum Lap Time: ',num2str(minTLap), '[s]']);
axis equal;
hold off;


xlabel('x(m)','fontweight','bold','fontsize',14)
ylabel('y(m)','fontweight','bold','fontsize',14)
title(sprintf('%s - Velocity Profile\nLap Time = %.2fs',name,time(end)),'fontsize',16)
