% Josh Sale
% Parse a CSV file from the RangeNet software
clear; close all; clc;
% Open the file 
fname = 'RangeNet_Static_3_1_19_002.csv';
fid = fopen(fname);

% Create an empty structure to hold all the data
uwbdata.t = [];
uwbdata.x = [];
uwbdata.y = [];
uwbdata.z = [];

% Read file line-by-line
str = fgetl(fid);
cnt = 1;
while (str ~= -1)
    % Split at commas
    C = strsplit(str,',');
    % Find the messages we are looking for
    if (length(C) == 23)
        if (strcmp(C{2},' LocEchoLastLocationExInfo') && ...
            ~strcmp(C{1},'Timestamp'))
            uwbdata.t(cnt) = str2double(C{1});
            uwbdata.x(cnt) = str2double(C{12});
            uwbdata.y(cnt) = str2double(C{13});
            uwbdata.z(cnt) = str2double(C{14});
            cnt = cnt+1;
        end
    end
    str = fgetl(fid);
end

% data.t=datetime(uwbdata.t, 'ConvertFrom', 'posixtime');

% Plot raw pulson x vs. y position
figure(1)
plot(uwbdata.x,uwbdata.y)
hold on
title('PulsON X vs. Y Position')
xlabel('X (mm)')
ylabel('Y (mm)')
hold off

uwbdata.t=uwbdata.t-8*3600; % Converts UWB time to UTC

%Plot pulson t vs. position
figure(2)
subplot(3,1,1)
plot(uwbdata.t,uwbdata.x)
hold on
title('PulsON T vs X Position')
xlabel('T')
ylabel('X')
% datetick('x',13,'keeplimits');
subplot(3,1,2)
plot(uwbdata.t,uwbdata.y)
title('PulsON T vs Y Position')
xlabel('T')
ylabel('Y')
% datetick('x',13,'keeplimits');
subplot(3,1,3)
plot(uwbdata.t,uwbdata.z)
title('PulsON T vs Z Position')
xlabel('T')
ylabel('Z')
% datetick('x',13,'keeplimits');
hold off

% Normalize data by removing initial vector [x0,y0,z0]
xdata=uwbdata.x-uwbdata.x(1);
ydata=uwbdata.y-uwbdata.y(1);
zdata=uwbdata.z-uwbdata.z(1);
t = datetime(uwbdata.t, 'ConvertFrom', 'posixtime')

%% Import and plot data run from VICON

load('data_03_01_19_2.mat')

vic_t=ans.Data(:,1); 
vic_x=ans.Data(:,2).*1000;
vic_y=ans.Data(:,3).*1000;
vic_z=-ans.Data(:,4).*1000;

vic_x1=vic_x-vic_x(1);
vic_y1=vic_y-vic_y(1);

mx=mean(vic_x);
my=mean(vic_y);
mz=mean(vic_z);
stdx=std(vic_x);
stdy=std(vic_y);
stdz=std(vic_z);

% Plot Vicon x vs. y position data
figure(3)
plot(vic_x,vic_y)
hold on
title('VICON 2-D Location')
hold on
xlabel('X (mm)')
ylabel('Y (mm)')

% Plot uwb & vicon t vs x data
figure(4)
plot(uwbdata.t,xdata,'b-')
hold on
plot(vic_t,vic_x1,'r-')
title('Vicon T vs. X')
xlabel('T')
ylabel('X')

%Plot pulson t vs. position
figure(5)
subplot(3,1,1)
plot(vic_t,vic_x)
hold on
title('VICON T vs X Position')
xlabel('T')
ylabel('X')
% datetick('x',13,'keeplimits');
subplot(3,1,2)
plot(vic_t,vic_y)
title('VICON T vs Y Position')
xlabel('T')
ylabel('Y')
% datetick('x',13,'keeplimits');
subplot(3,1,3)
plot(vic_t,vic_z)
title('VICON T vs Z Position')
xlabel('T')
ylabel('Z')
% datetick('x',13,'keeplimits');
hold off

% Plot Vicon & UWB x vs. y position data on same axes
figure(6)
plot(vic_x,vic_y,'r-')
hold on
title('2-D Location')
plot(uwbdata.x,uwbdata.y,'b-')
hold on
legend('Vicon Data', 'PulsON Data','Location','southeast');
xlabel('X (mm)')
ylabel('Y (mm)')
hold off

% Translate and rotate vicon data to overlay on UWB data
psi=tan((xdata(1000)-xdata(1))/ydata(1000)-ydata(1));
R=[cos(psi) -sin(psi);sin(psi) cos(psi)];
vic_rot=(R*[vic_x1'; vic_y1'])';
vic_rotx=vic_rot(:,1);
vic_roty=vic_rot(:,2);

% Plot vicon & UWB x vs. y position data on same axes
figure(7)
plot(vic_rotx,vic_roty,'r-')
hold on
title('2-D Location')
plot(xdata,ydata,'b-')
legend('Vicon Data', 'PulsON Data','Location','southeast');
xlabel('X (mm)')
ylabel('Y (mm)')
