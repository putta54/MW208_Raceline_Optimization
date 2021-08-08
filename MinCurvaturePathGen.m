%% Processing  track data

%track data - first point repeated
data = track_silverstone;

% x,y and track width data 
x =  data(:,1);
y =  data(:,2);
twr = data(:,3);
twl = data(:,4);

% interpolate data to get finer curve
t = [0; cumsum(hypot(diff(x),diff(y)))];
t1 = linspace(0,t(end),1500);
xt = spline(t,x,t1)';
yt = spline(t,y,t1)';
twrt = spline(t,twr,t1)';
twlt = spline(t,twl,t1)';
% xt = x;
% yt = y;

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

% % plot segments
% figure
% hold on
% for i=1:numel(xout)
%     plot([xin(i) xout(i)], [yin(i) yout(i)])
% end

% form delta matrices
delx = xout - xin;
dely = yout - yin;

%% Matrix Definition 

% number of segments
n = numel(delx);

% preallocation
H = zeros(n);
B = zeros(size(delx)).';

% formation of H matrix (nxn)
for i=2:n-1
    
    % first row
    H(i-1,i-1) = H(i-1,i-1) + delx(i-1)^2         + dely(i-1)^2;
    H(i-1,i)   = H(i-1,i)   - 2*delx(i-1)*delx(i) - 2*dely(i-1)*dely(i);
    H(i-1,i+1) = H(i-1,i+1) + delx(i-1)*delx(i+1) + dely(i-1)*dely(i+1);
    
    %second row
    H(i,i-1)   = H(i,i-1)   - 2*delx(i-1)*delx(i) - 2*dely(i-1)*dely(i);
    H(i,i)     = H(i,i )    + 4*delx(i)^2         + 4*dely(i)^2;
    H(i,i+1)   = H(i,i+1)   - 2*delx(i)*delx(i+1) - 2*dely(i)*dely(i+1);
    
    % third row
    H(i+1,i-1) = H(i+1,i-1) + delx(i-1)*delx(i+1) + dely(i-1)*dely(i+1);
    H(i+1,i)   = H(i+1,i)   - 2*delx(i)*delx(i+1) - 2*dely(i)*dely(i+1);
    H(i+1,i+1) = H(i+1,i+1) + delx(i+1)^2         + dely(i+1)^2;
    
end

% formation of B matrix (1xn)
for i=2:n-1
    
    B(1,i-1) = B(1,i-1) + 2*(xin(i+1)+xin(i-1)-2*xin(i))*delx(i-1) + 2*(yin(i+1)+yin(i-1)-2*yin(i))*dely(i-1);
    B(1,i)   = B(1,i)   - 4*(xin(i+1)+xin(i-1)-2*xin(i))*delx(i)   - 4*(yin(i+1)+yin(i-1)-2*yin(i))*dely(i);
    B(1,i+1) = B(1,i+1) + 2*(xin(i+1)+xin(i-1)-2*xin(i))*delx(i+1) + 2*(yin(i+1)+yin(i-1)-2*yin(i))*dely(i+1);
    
end

% define constraints
lb = zeros(n,1);
ub = ones(size(lb));

% if start and end points are the same
Aeq      =   zeros(1,n);
Aeq(1)   =   1;
Aeq(end) =   -1;
beq      =   0;
    
%% Solver

options = optimoptions('quadprog','Display','iter');
[resMCP,fval,exitflag,output] = quadprog(2*H,B',[],[],Aeq,beq,lb,ub,[],options);

%% Plotting results

% co-ordinates for the resultant curve
xresMCP = zeros(size(xt));
yresMCP = zeros(size(xt));

for i = 1:numel(xt)
    xresMCP(i) = xin(i)+resMCP(i)*delx(i);
    yresMCP(i) = yin(i)+resMCP(i)*dely(i);
end

figure
plot(xresMCP,yresMCP,'color','r','linew',2)
hold on

% plot starting line
plot([xin(1) xout(1)], [yin(1) yout(1)],'color','b','linew',2)
% plot([xin(2) xout(2)], [yin(2) yout(2)],'color','k','linew',2)

% plot reference line
plot(x,y,'--')
hold on

% plot inner track
plot(xin,yin,'color','k')

%plot outer track
plot(xout,yout,'color','k')
hold off
axis equal

xlabel('x(m)','fontweight','bold','fontsize',14)
ylabel('y(m)','fontweight','bold','fontsize',14)
title('Silverstone, UK (F1) - Minimum Curvature Trajectory','fontsize',16)