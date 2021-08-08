%% Input data form

% track data form = [x y track_width_to_the_right(+ve) track_width_to_the_left(+ve)]
% name = 'name_of_track'

function [xresSP, yresSP] = func_SP(track,name)
%% Processing  track data

% track data - first point repeated
data = track;

% x, y and track width data 
x =  data(:,1);
y =  data(:,2);
twr = data(:,3);
twl = data(:,4);

% interpolate data to get finer curve with equal distances between each segment
pathXY = [x y];
stepLengths = sqrt(sum(diff(pathXY,[],1).^2,2));
stepLengths = [0; stepLengths]; % add the starting point
cumulativeLen = cumsum(stepLengths);
finalStepLocs = linspace(0,cumulativeLen(end), 1500);
finalPathXY = interp1(cumulativeLen, pathXY, finalStepLocs);
xt = finalPathXY(:,1);
yt = finalPathXY(:,2);
twrt = interp1(cumulativeLen, twr, finalStepLocs,'spline')';
twlt = interp1(cumulativeLen, twl, finalStepLocs,'spline')';

% normal direction for each vertex
dx = gradient(xt);
dy = gradient(yt);
dL = hypot(dx,dy);

% offset curve - anonymous function
xoff = @(a) -a*dy./dL + xt;
yoff = @(a)  a*dx./dL + yt;

% plot reference line
plot(xt,yt,'g')
hold on

% offset data
offset = [-twrt twlt];
for i = 1:numel(xt)
    xin = xoff(offset(i,1));      % get inner offset curve
    yin = yoff(offset(i,1));
    
    xout  = xoff(offset(i,2));      % get outer offset curve
    yout  = yoff(offset(i,2));
end

% plot inner track
plot(xin,yin,'color','b','linew',2)

% plot outer track
plot(xout,yout,'color','r','linew',2)
hold off
axis equal

xlabel('x(m)','fontweight','bold','fontsize',14)
ylabel('y(m)','fontweight','bold','fontsize',14)
title(sprintf(name),'fontsize',16)

% % plot segments
% figure
% hold on
% for i=1:numel(xout)
%     plot([xin(i) xout(i)], [yin(i) yout(i)])
% end

% form delta matrices
delx = xout - xin;
dely = yout - yin;

%% Matrix Definition (H and B)

% number of segments
n = numel(delx);

% preallocation
H = zeros(n);
B = zeros(size(delx)).';

% formation of H matrix (nxn)
for i=1:n-1
    
    H(i,i)     = H(i,i)     + delx(i)^2          + dely(i)^2;
    H(i+1,i)   = H(i+1,i)   - delx(i)*delx(i+1)  - dely(i)*dely(i+1);
    H(i,i+1)   = H(i,i+1)   - delx(i)*delx(i+1)  - dely(i)*dely(i+1);
    H(i+1,i+1) = H(i+1,i+1) + delx(i+1)^2        + dely(i+1)^2;
    
end

% formation of B matrix (1xn)
for i=1:n-1
    B(1,i)   = B(1,i)   - 2*(xin(i+1)-xin(i))*delx(i)   - 2*(yin(i+1)-yin(i))*dely(i);
    B(1,i+1) = B(1,i+1) + 2*(xin(i+1)-xin(i))*delx(i+1) + 2*(yin(i+1)-yin(i))*dely(i+1);
end

% define boundary constraints
lb = zeros(n,1);
ub = ones(size(lb));

% forcing start and end points to be the same
Aeq      =   zeros(1,n);
Aeq(1)   =   1;
Aeq(end) =   -1;
beq      =   0;
    
%% QP Solver

options = optimoptions('quadprog','Display','iter');
[resSP,fval,exitflag,output] = quadprog(2*H,B,[],[],Aeq,beq,lb,ub,[],options);

%% Plotting results

% co-ordinates for the resultant curve
xresSP = zeros(size(xt));
yresSP = zeros(size(xt));

for i = 1:numel(xt)
    xresSP(i) = xin(i)+resSP(i)*delx(i);
    yresSP(i) = yin(i)+resSP(i)*dely(i);
end

% plot shortest path
figure
plot(xresSP,yresSP,'color','r','linew',2)
hold on

% plot reference line
plot(x,y,'--')
hold on

% plot inner track
plot(xin,yin,'color','k')

% plot outer track
plot(xout,yout,'color','k')
hold off
axis equal

xlabel('x(m)','fontweight','bold','fontsize',14)
ylabel('y(m)','fontweight','bold','fontsize',14)
title(sprintf(name,'- Shortest Path Trajectory'),'fontsize',16)
