% crazyflie load usd card log csv file

% read file
csvfilename = '2024-04-04_log07.csv';
T = readtable(csvfilename);

% get data from table
time = table2array(T(:,1))'*1e-3;
pos = table2array(T(:,2:4))';
vel = table2array(T(:,5:7))';
% acc = table2array(T(:,8:10))';
lbd = table2array(T(:,8:10))'*pi/180;
om = table2array(T(:,11:13))'*pi/180;
pos_ref = table2array(T(:,14:16))';
yaw_ref = table2array(T(:,17))';
motors = table2array(T(:,18:21))';

Ts = mean(diff(time));

px1_id = pos(1,9348:9601)';
py1_id = pos(2,9348:9601)';
pz1_id = pos(3,9348:9601)';

fi1_id = lbd(1,9348:9601)';
theta1_id = lbd(2,9348:9601)';
psi1_id = lbd(3,9348:9601)';

px2_id = pos(1,19108:19361)';
py2_id = pos(2,19108:19361)';
pz2_id = pos(3,19108:19361)';

fi2_id = lbd(1,19108:19361)';
theta2_id = lbd(2,19108:19361)';
psi2_id = lbd(3,19108:19361)';

% convert date to print format
t = time - time(1);
x = [pos;vel;lbd;om];
x_ref = [pos_ref;0*vel;lbd*0;om*0];
x_ref(9,:) = yaw_ref;
u = motors/65535.0;

T1_id = u(1,9348:9601)'; 
T2_id = u(1,19108:19361)'; 

% plot data
initPlots;
vehicle3d_ref_show_data(t,x,u,x_ref);