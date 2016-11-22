clear all;

% xOrig = [1,1.5,3,7,8,9,10,11,12,19,20,21,22,23,24,27,28,29,30,31,32,33,34,35,36,37];
% yOrig = [1,1.5,3,3,2,1,0,-1,-2,3,4,4.5, 5, 5, 5, 2, 2, 2, 2, 2.5, 4, 1, 3, 2, 3, 3];

%xOrig = [0.5,1,2,3,4,8,9,10,11];
%yOrig = [0.5,1,3,5,7,33,30,27,24];

%xOrig = [1,4,5,6,7,9,10,11,12,14];
%yOrig = [1,2,3,4,5,6,7,8,9,10];

%xOrig = [1,2,3,4,5,11,12,13,14,15];
%yOrig = [1,1,1,1,1,3,3,3,3,3];

xOrig = [1,2,3,4,5,9,10,11,12,13,14,15,21,22,23,24];
yOrig = [1,2,8,7,6,1,1,2,2,3,4,5,2,-3,-10,10];

%add noise to data
for i = 1:size(xOrig,2)
    yOrig(i) = yOrig(i) + rand()*0.8;
end

x = xOrig;
y = yOrig;

fprintf('size of xOrig=%d \n', size(xOrig,2));
fprintf('size of yOrig=%d \n' , size(yOrig,2));

%find shortest period
shortestPeriod = 9999;
for i = 1:(size(x,2)-1)
    diff = x(i+1)-x(i);
    if diff < shortestPeriod
        shortestPeriod = diff;
    end
end

%New trajectory
xn = [];
yn = [];

xn = [xn, x(1)];

speedLimit = 2;
timeCompensation = 0;
stretchFactor = 1.2;

%find the sections that exceed the joint speed limit
for i = 1:(size(x,2)-1)
    diff = x(i+1)-x(i);
    vel = abs( (y(i+1)-y(i)) / diff );
    
    if vel >= speedLimit
        fprintf('speed exceeded at i=%d, speed=%f\n', i, vel)
        newDiff = (vel/speedLimit)*diff*stretchFactor;
        timeCompensation = timeCompensation + (newDiff - diff);
        fprintf('timeCompensation=%f\n', timeCompensation)
    end
    
    xn = [xn, x(i+1) + timeCompensation];    
end

x = xn;
xn = [];

%start the interpolation
disp('starting interpolation search ')

%start looking for larger gaps first and then reduce the gap threshold per
%eac cycle
for k = linspace(5,2,3)
    for i = 1:(size(x,2)-1)
        %fprintf('  i=%d \n', i)
        diff = x(i+1)-x(i);
        %fprintf('  diff=%f \n', diff)
        
        % 
        if diff > k*shortestPeriod
            %disp('    found a gap larger than shortestPeriod')
            %Get the velocities
            vel1 = 0;
            vel2 = 0;

            %If first or last period, then have same velocities
            if (i == 1) | (i == (size(x,2)-1))
                disp('    first/last period');
                vel1 = (y(i+1)-y(i)) / diff;
                vel2 = vel1;       
            else
                disp('    intermeidate period');
                vel1 = (y(i)-y(i-1)) / (x(i)-x(i-1));
                vel2 = (y(i+2)-y(i+1)) / (x(i+2)-x(i+1));
            end

            fprintf('    vel1=%f \n', vel1);
            fprintf('    vel2=%f \n\n', vel2);

            %initialize the variables
            acceleration = (vel2 - vel1) / diff;
            stepsFloat = diff / shortestPeriod;
            stepsInt = floor(stepsFloat);
            dt = diff / stepsInt;

            %first, add the starting point
            xn = [xn, x(i)];
            yn = [yn, y(i)];

            %start filling the area in between
            for j = 1:(stepsInt - 1)
                %fprintf('      j=%d \n', j)
                new_x = x(i) + dt*j;
                traj1 = y(i) + vel1*dt*j + (acceleration/2)*((dt*j)^2);
                traj2 = y(i+1) - vel2*dt*(stepsInt - j) + (acceleration/2)*((dt*(stepsInt - j))^2);

                %blend = (dt*j/diff)^1;
                blend = (1 + cos(pi*dt*j/diff))/2;
                new_y = blend*traj1 + (1 - blend)*traj2;
                %new_y = traj2;
                
                %fprintf('      new_x=%f \n', new_x);
                %fprintf('      new_y=%f \n', new_y);

                xn = [xn, new_x];
                yn = [yn, new_y];
            end

        else
            xn = [xn, x(i)];
            yn = [yn, y(i)];
        end
    end
    
    xn = [xn, x(end)];
    yn = [yn, y(end)];
    
    x = xn;
    y = yn;
    
    xn = [];
    yn = [];
end

% Gaussian filtering
gaussFilter = [0.06136,	0.24477, 0.38774, 0.24477, 0.06136];
kernelcheck = sum(gaussFilter)
startIndex = floor(size(gaussFilter,2)/2) + 1

xg = x;
yg = [y(1), y(2)];

for i = 3:(size(x,2)-2)
    fprintf('i=%d \n', i);
    filtered = 0;
    for j = 1:size(gaussFilter,2)
        filtered = filtered + gaussFilter(j)*y(i - startIndex + j);
    end
    yg = [yg, filtered];
end

yg = [yg, y(size(x,2) - 1), y(size(x,2))];

fprintf('size of xOrig=%d \n', size(xOrig,2));
fprintf('size of yOrig=%d \n' , size(yOrig,2));

fprintf('size of x=%d \n', size(x,2));
fprintf('size of y=%d \n' , size(y,2));

fprintf('size of xg=%d \n', size(xg,2));
fprintf('size of yg=%d \n' , size(yg,2));

plot(x, y, '-xb')
hold on;
plot(xOrig, yOrig, '-or')
hold on;
plot(xg, yg, '-+')

legend('interpolated','original','interp + gauss filt.','Location','southwest')

 
    
    