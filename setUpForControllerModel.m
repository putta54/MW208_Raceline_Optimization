
%% define reference points
xRef = trajMCP(:,2);
yRef = trajMCP(:,1);
refPose = [xRef yRef];

%% define vehicle parameters used in the models
X_o = xRef(1); % initial vehicle position in x direction
Y_o = yRef(1); % initial vehicle position in y direction

%% calculating reference pose vectors
% Based on how far the vehicle travels, the pose is generated using 1-D
% lookup tables.

% calculate distance vector
distancematrix = squareform(pdist(refPose));
distancesteps = zeros(length(refPose)-1,1);
for i = 2:length(refPose)
    distancesteps(i-1,1) = distancematrix(i,i-1);
end
totalDistance = sum(distancesteps); % Total traveled distance
distbp = cumsum([0; distancesteps]); % Distance for each waypoint
gradbp = linspace(0,totalDistance,2000); % Linearize distance

% linearize X and Y vectors based on distance
xRef2 = interp1(distbp,xRef,gradbp);
yRef2 = interp1(distbp,yRef,gradbp);
yRef2s = smooth(gradbp,yRef2);
xRef2s = smooth(gradbp,xRef2);

%% calculate curvature vector
curvature = getCurvature(xRef2,yRef2);

%% calculate theta vector
% theta = orientation angle of the path at reference points 
thetaRef = zeros(length(gradbp),1);
for i = 2:length(gradbp)
    thetaRef(i,1) = atan2d((yRef2(i)-yRef2(i-1)),(xRef2(i)-xRef2(i-1)));
end
thetaRefs = smooth(gradbp,thetaRef); % smoothing of theta


%% Removing abrupt changes in theta 

% Please refer to the "README.html" file understand this section with images. 

% The built-in function, "findchangepts" finds abrupt changes in the signal.
% When using atan2, if theta is tending towards |180| degrees, the interpolation might
% generate additional points leading to a fluctuating signal between -180
% degrees and 180 degrees. Hence, we follow these steps to smoothen the
% signal:
% -> Find abrupt changes in the signal using findchangepts
% -> Visualize the plot and find the change points where the change is
% causing significant deviation
% -> Define the regions where the signal has to be replaced by a
% neighboring smooth value
% Note: It is a manual process. So it's recommended to visualize the data and tune the
% changepoints and regions to remove the abrupt changes from the signal
% 
% If there is no fluctuation in theta i.e. if -180<theta<=180, consider removing this section and
% use thetaRefs as input to the Simulink model lookup table block named, "Theta Ref"

idx = 10; % number of change points to be visaulized
[ipt, residual] = findchangepts(thetaRefs, 'Statistic','linear','MaxNumChanges',idx);

%% Uncomment this section to visualize the change points
figure % figure to visualize the change points
findchangepts(thetaRefs, 'Statistic','linear','MaxNumChanges',idx)

% select the changepoints 
gradbp0 = ipt(1);
gradbp1 = ipt(2);
gradbp2 = max(ipt);
thetaRefab = zeros(length(gradbp),1); % theta after removing the abrupt changes 
% define the regions and remove the abrupt changes from the signal
for i = 1:length(gradbp)
    if (gradbp(i) == gradbp(gradbp1))
        thetaRefab(i) = thetaRefs(gradbp1);
    else
        thetaRefab(i) = thetaRefs(i);
    end
end
psi_o = thetaRefab(1)*(pi/180); % initial yaw angle

%% plot to visualize the difference between thetaRefs and thetaRefab
figure
plot(gradbp,thetaRefab,'Linewidth',5,'Color','r')
hold on
plot(gradbp,thetaRefs,'b')
xlabel('distance (m)')
ylabel('theta (deg)')

%% create direction vector
direction = ones(length(gradbp),1);

%% Curvature Function

function curvature = getCurvature(xRef,yRef)
% Calculate gradient by the gradient of the X and Y vectors
DX = gradient(xRef);
D2X = gradient(DX);
DY = gradient(yRef);
D2Y = gradient(DY);
curvature = (DX.*D2Y - DY.*D2X) ./(DX.^2+DY.^2).^(3/2);
end
