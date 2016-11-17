clear all;

xOrig = [1,1.5,3,7,8,9,10,11,12,19,20,21,22,23,24,27,28,29,30];
yOrig = [1,1.5,3,3,2,1,0,-1,-2,3,4,4.5, 5, 5, 5, 2, 2, 2, 2];

%xOrig = [0.5,1,2,3,4,8,9,10,11];
%yOrig = [0.5,1,3,5,7,33,30,27,24];

%xOrig = [1,4,5,6,7,9,10,11,12,14];
%yOrig = [1,2,3,4,5,6,7,8,9,10];

%xOrig = [1,2,3,4,5,11,12,13,14,15];
%yOrig = [1,1,1,1,1,3,3,3,3,3];


%add noise to data
for i = 1:size(xOrig,2)
    yOrig(i) = yOrig(i) + rand()*1.5;
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

%start the interpolation
disp('starting interpolation search ')

%start looking for larger gaps first and then reduce the gap threshold per
%eac cycle
for k = linspace(5,1,3)
    for i = 1:(size(x,2)-1)
        fprintf('  i=%d \n', i)
        diff = x(i+1)-x(i);
        fprintf('  diff=%f \n', diff)
        
        % 
        if diff > k*shortestPeriod
            disp('    found a gap larger than shortestPeriod')
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
                fprintf('      j=%d \n', j)
                new_x = x(i) + dt*j;
                traj1 = y(i) + vel1*dt*j + (acceleration/2)*((dt*j)^2);
                traj2 = y(i+1) - vel2*dt*(stepsInt - j) + (acceleration/2)*((dt*(stepsInt - j))^2);

                %blend = (dt*j/diff)^1;
                blend = (1 + cos(pi*dt*j/diff))/2;
                new_y = blend*traj1 + (1 - blend)*traj2;
                %new_y = traj2;
                
                fprintf('      new_x=%f \n', new_x);
                fprintf('      new_y=%f \n', new_y);

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
    
fprintf('size of xOrig=%d \n', size(xOrig,2));
fprintf('size of yOrig=%d \n' , size(yOrig,2));

fprintf('size of xn=%d \n', size(xn,2));
fprintf('size of yn=%d \n' , size(yn,2));

plot(x, y, '-x')
hold
plot(xOrig, yOrig, '-o')


 
    
    