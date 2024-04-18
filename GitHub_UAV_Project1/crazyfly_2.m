clc
clear all;
initPlots;

% Simulation parameters
nx = 12 ;
ny = 12 ;
x0 = [0; 0; 0; 0 ; 0 ; 0 ; 0; 0; 0; 0; 0; 0] ;
Dt = 0.001 ;
t = 0:Dt:5 ;
g = 9.81 ;
rho = 1.225 ;
Tmax = .7848/4 ;
cq =.10 ;
ct = 0.15 ;
lx = 0.65 ;
ly = 0.065 ;
lz = 0.029 ;
l = 0.046 ;
Ax = ly * lz ;
Ay = lx * lz ;
Az = lx * ly ;
m = 0.028 ;
cd = 2.3 ;
betax = .5*cd*rho*Ax ;
betay = .5*cd*rho*Ay ;
betaz = .5*cd*rho*Az ;
D = diag([ betax , betay , betaz ]) ;
zI = [0 ; 0 ; -1] ;
ang = 45*pi/180 ;
Jx = 1/12* m*(ly^2+lz^2) ;
Jy = 1/12* m*(lx^2+lz^2) ;
Jz = 1/12* m*(lx^2+ly^2) ;
J = diag([ Jx , Jy , Jz ]) ;

p1 = l*[ cos(ang) ; -sin(ang) ; 0 ]  ;
p2 = l*[ -cos(ang) ; -sin(ang) ; 0 ]  ;
p3 = l*[ -cos(ang) ; sin(ang) ; 0 ]  ;
p4 = l*[ cos(ang) ; sin(ang) ; 0 ]   ;

T1 = 0.4*Tmax ;
T2 = T1 ;
T3 = T1 ;
T4 = T1 ;

fp1 = T1*[0 ; 0 ; 1] ;
fp2 = T2*[0 ; 0 ; 1] ;
fp3 = T3*[0 ; 0 ; 1] ;
fp4 = T4*[0 ; 0 ; 1];

syms Q1 Q2  Q3 Q4
n1 = Q1*[0 ; 0 ; -1] ;
n2 = Q2*[0 ; 0 ; 1] ;
n3 = Q3*[0 ; 0 ; -1] ;
n4 = Q4*[0 ; 0 ; 1] ;

np1 = n1 + skew(p1)*fp1 ;
np2 = n2 + skew(p2)*fp2 ;
np3 = n3 + skew(p3)*fp3 ;
np4 = n4 + skew(p4)*fp4 ;

fp = fp1 + fp2 + fp3 + fp4 ;

np = np1 + np2 + np3 + np4 ;

a = sqrt(2)/2 * l ;

T1 = 0.4*Tmax ;
T2 = T1 ;
T3 = T1 ;
T4 = T1 ;

u_NL = [ 1 , 1 , 1 , 1 ; 
        -a , -a , a , a  ;
        -a, a , -a , a ;  
        -cq/ct , cq/ct , -cq/ct , cq/ct ] *[T1 ; T2 ; T3 ; T4]* ones(size(t));


% simulate nonlinear system
Nsim = length(t);
x = zeros(nx,Nsim);
y = zeros(ny,Nsim);
x(:,1) = x0;

for k = 1:Nsim
    % prepare variables:
    p   = x(1:3,k);
    v = x(4:6,k);
    lbd  = x(7:9,k);
    om  = x(10:12,k);
    R = Euler2R(lbd);
    T = u_NL(1,k);
    np = u_NL(2:4,k);
    
    % compute state derivative:
    p_dot = Euler2R(lbd)*v;
    lbd_dot = Euler2Q(lbd)*om;
    v_dot = -skew(om)*v + g*R'*zI - (1/m)*R*D*R'* v .* abs(v) + (1/m)*fp ; 
    om_dot = -inv(J)*skew(om)*J*om + inv(J)*np;
    x_dot = [p_dot;v_dot;lbd_dot;om_dot];

    % integrate state
    x(:,k+1) = x(:,k) + Dt*x_dot;

    % compute current output:
    y(:,k) = x(:,k);
end


figure(1);
plot(t,y(3,:), '-.',t,y(6,:), '-.');
grid on;
xlabel('Time [s]');
legend( '$$p_z$$','$$v_z$$');
title('Nonlinear Model Simulation');

%Operation Point 2
px = 0 ;
py = 0 ;
pz =  3 ; 
fi = 0 ;
theta = .5*pi/180 ;
psi = 0 ;
vy = 0 ;
vz = 0 ;
wx =0 ;
wy =0 ;
wz = 0 ;


syms T vx
eqn1 = 4*T*cos(theta) - m*g == 0  ;   %Fz
T = double (solve(eqn1,T))


eqn2 = -betax*vx*abs(vx) +4*T*sin(theta) == 0 ; %Fy
vx = double (solve(eqn2,vx))




T1 = T ;
T2 = T1 ;
T3 = T1 ;
T4 = T1 ;

xop = [px ; py ; pz ; vx ; vy ; vz ; fi ; theta ; psi ; wx ; wy ; wz ] ;

A1 = zeros (3,3) ; 
%p/p

A2 = [ vy*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)) + vz*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)), vz*cos(fi)*cos(psi)*cos(theta) - vx*cos(psi)*sin(theta) + vy*cos(psi)*cos(theta)*sin(fi), vz*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)) - vy*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)) - vx*cos(theta)*sin(psi);
      - vy*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)) -  vz*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)), vz*cos(fi)*cos(theta)*sin(psi) - vx*sin(psi)*sin(theta) + vy*cos(theta)*sin(fi)*sin(psi), vz*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)) - vy*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)) + vx*cos(psi)*cos(theta);
        vy*cos(fi)*cos(theta) - vz*cos(theta)*sin(fi), - vx*cos(theta) - vz*cos(fi)*sin(theta) - vy*sin(fi)*sin(theta),0] ;
%p/lbd 
    
A3 = [ cos(theta)*cos(psi) , sin(fi)*sin(theta)*cos(psi)-cos(fi)*sin(psi) , cos(fi)*sin(theta)*cos(psi)+sin(fi)*sin(psi) ;
       cos(theta)*sin(psi)  , sin(fi)*sin(theta)*sin(psi)+cos(fi)*cos(psi) , cos(fi)*sin(theta)*sin(psi)-sin(fi)*cos(psi) ;
       -sin(theta)          , sin(fi)*cos(theta)                           , cos(fi)*cos(theta) ] ; 
%p/v 

A4 = zeros (3,3) ;
%p/om

A5 = zeros (3,3) ;
%lbd/p
 
A6 = [ wy*cos(fi)*tan(theta) - wz*sin(fi)*tan(theta)    , wz*cos(fi)*(tan(theta)^2 + 1) + wy*sin(fi)*(tan(theta)^2 + 1)               , 0 ;                                          
        - wz*cos(fi) - wy*sin(fi)                                      , 0                                                              , 0 ;
       (wy*cos(fi))/cos(theta) - (wz*sin(fi))/cos(theta) ,(wz*cos(fi)*sin(theta))/cos(theta)^2 + (wy*sin(fi)*sin(theta))/cos(theta)^2  , 0 ] ;                         
%lbd/lbd

A7 = zeros (3,3) ;
%lbd/v
 
A8 = [ 1 , sin(fi)*tan(theta) , cos(fi)*tan(theta) ;
       0  , cos(fi)            , -sin(fi) ;
       0  , sin(fi)/cos(theta) , cos(fi)/cos(theta) ] ;
%lbd/om

A9 = zeros (3,3) ;
%v/p

A10 = [ abs(vx)*(vx*((betay*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta)))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m + (betay*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta)))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m - (betaz*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta)))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m - (betaz*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta)))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m) - vy*((betay*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m - (betaz*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m + (betay*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m - (betaz*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m) + vz*((betay*cos(conj(fi))*cos(conj(theta))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m - (betaz*cos(conj(fi))*cos(conj(theta))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m - (betay*cos(conj(theta))*sin(conj(fi))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m + (betaz*cos(conj(theta))*sin(conj(fi))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m)) ,  g*cos(conj(theta)) + abs(vx)*(vx*((betay*cos(psi)*cos(theta)*sin(fi)*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta))))/m - (betaz*cos(fi)*cos(psi)*cos(theta)*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta))))/m - (betaz*cos(conj(fi))*cos(conj(psi))*cos(conj(theta))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m + (betay*cos(conj(psi))*cos(conj(theta))*sin(conj(fi))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m + (betax*cos(conj(psi))*cos(conj(theta))*cos(psi)*sin(theta))/m + (betax*cos(conj(psi))*sin(conj(theta))*cos(psi)*cos(theta))/m) - vz*((betay*sin(conj(fi))*sin(conj(theta))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m - (betax*cos(conj(theta))*cos(psi)*cos(theta))/m - (betaz*cos(conj(fi))*sin(conj(theta))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m + (betax*sin(conj(theta))*cos(psi)*sin(theta))/m + (betaz*cos(conj(fi))*cos(conj(theta))*cos(fi)*cos(psi)*cos(theta))/m + (betay*cos(conj(theta))*sin(conj(fi))*cos(psi)*cos(theta)*sin(fi))/m) + vy*((betaz*cos(fi)*cos(psi)*cos(theta)*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m - (betay*cos(psi)*cos(theta)*sin(fi)*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m - (betaz*cos(conj(fi))*cos(conj(theta))*sin(conj(psi))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m + (betay*cos(conj(theta))*sin(conj(fi))*sin(conj(psi))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m + (betax*cos(conj(theta))*sin(conj(psi))*cos(psi)*sin(theta))/m + (betax*sin(conj(psi))*sin(conj(theta))*cos(psi)*cos(theta))/m))  , -abs(vx)*(vz*((betaz*cos(conj(fi))*cos(conj(theta))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m - (betay*cos(conj(theta))*sin(conj(fi))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m + (betax*sin(conj(theta))*cos(theta)*sin(psi))/m) + vx*((betay*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m + (betay*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta)))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m + (betaz*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta)))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m + (betaz*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m - (betax*cos(conj(psi))*cos(conj(theta))*cos(theta)*sin(psi))/m - (betax*cos(conj(theta))*sin(conj(psi))*cos(psi)*cos(theta))/m) - vy*((betay*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m - (betaz*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta)))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m - (betay*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta)))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m + (betaz*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m - (betax*cos(conj(psi))*cos(conj(theta))*cos(psi)*cos(theta))/m + (betax*cos(conj(theta))*sin(conj(psi))*cos(theta)*sin(psi))/m))                  ;
       - abs(vy)*(vx*((betay*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta)))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m - (betaz*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta)))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m + (betay*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta)))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m - (betaz*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta)))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m) - vy*((betay*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m + (betay*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m - (betaz*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m - (betaz*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m) + vz*((betay*cos(conj(fi))*cos(conj(theta))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m - (betaz*cos(conj(fi))*cos(conj(theta))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m - (betay*cos(conj(theta))*sin(conj(fi))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m + (betaz*cos(conj(theta))*sin(conj(fi))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m)) - g*cos(conj(fi))*cos(conj(theta)) , abs(vy)*(vy*((betaz*cos(fi)*cos(theta)*sin(psi)*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m - (betay*cos(theta)*sin(fi)*sin(psi)*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m + (betaz*cos(conj(fi))*cos(conj(theta))*sin(conj(psi))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m - (betay*cos(conj(theta))*sin(conj(fi))*sin(conj(psi))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m + (betax*cos(conj(theta))*sin(conj(psi))*sin(psi)*sin(theta))/m + (betax*sin(conj(psi))*sin(conj(theta))*cos(theta)*sin(psi))/m) - vz*((betaz*cos(conj(fi))*sin(conj(theta))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m - (betay*sin(conj(fi))*sin(conj(theta))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m - (betax*cos(conj(theta))*cos(theta)*sin(psi))/m + (betax*sin(conj(theta))*sin(psi)*sin(theta))/m + (betaz*cos(conj(fi))*cos(conj(theta))*cos(fi)*cos(theta)*sin(psi))/m + (betay*cos(conj(theta))*sin(conj(fi))*cos(theta)*sin(fi)*sin(psi))/m) + vx*((betay*cos(theta)*sin(fi)*sin(psi)*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta))))/m - (betaz*cos(fi)*cos(theta)*sin(psi)*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta))))/m + (betaz*cos(conj(fi))*cos(conj(psi))*cos(conj(theta))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m - (betay*cos(conj(psi))*cos(conj(theta))*sin(conj(fi))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m + (betax*cos(conj(psi))*cos(conj(theta))*sin(psi)*sin(theta))/m + (betax*cos(conj(psi))*sin(conj(theta))*cos(theta)*sin(psi))/m)) + g*sin(conj(fi))*sin(conj(theta))  ,  abs(vy)*(vz*((betay*cos(conj(theta))*sin(conj(fi))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m - (betaz*cos(conj(fi))*cos(conj(theta))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m + (betax*sin(conj(theta))*cos(psi)*cos(theta))/m) + vx*((betay*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m - (betaz*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta)))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m - (betay*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta)))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m + (betaz*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m - (betax*cos(conj(psi))*cos(conj(theta))*cos(psi)*cos(theta))/m + (betax*cos(conj(theta))*sin(conj(psi))*cos(theta)*sin(psi))/m) + vy*((betay*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m + (betay*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta)))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m + (betaz*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta)))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m + (betaz*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m - (betax*cos(conj(psi))*cos(conj(theta))*cos(theta)*sin(psi))/m - (betax*cos(conj(theta))*sin(conj(psi))*cos(psi)*cos(theta))/m))                  ;
       g*cos(conj(theta))*sin(conj(fi)) - abs(vz)*(vy*((betay*cos(fi)*cos(theta)*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m - (betaz*cos(fi)*cos(theta)*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m - (betay*cos(theta)*sin(fi)*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m + (betaz*cos(theta)*sin(fi)*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m) - vx*((betay*cos(fi)*cos(theta)*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta))))/m - (betaz*cos(fi)*cos(theta)*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta))))/m - (betay*cos(theta)*sin(fi)*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta))))/m + (betaz*cos(theta)*sin(fi)*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta))))/m) + vz*((betay*cos(conj(fi))*cos(conj(theta))*cos(theta)*sin(fi))/m + (betay*cos(conj(theta))*sin(conj(fi))*cos(fi)*cos(theta))/m - (betaz*cos(conj(fi))*cos(conj(theta))*cos(theta)*sin(fi))/m - (betaz*cos(conj(theta))*sin(conj(fi))*cos(fi)*cos(theta))/m)) , g*cos(conj(fi))*sin(conj(theta)) - abs(vz)*(vy*((betax*sin(conj(psi))*sin(conj(theta))*sin(theta))/m + (betaz*cos(fi)*sin(theta)*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m - (betay*sin(fi)*sin(theta)*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m - (betax*cos(conj(theta))*sin(conj(psi))*cos(theta))/m + (betaz*cos(conj(fi))*cos(conj(theta))*sin(conj(psi))*cos(fi)*cos(theta))/m + (betay*cos(conj(theta))*sin(conj(fi))*sin(conj(psi))*cos(theta)*sin(fi))/m) - vz*((betaz*cos(conj(fi))*cos(conj(theta))*cos(fi)*sin(theta))/m - (betax*sin(conj(theta))*cos(theta))/m - (betax*cos(conj(theta))*sin(theta))/m + (betaz*cos(conj(fi))*sin(conj(theta))*cos(fi)*cos(theta))/m + (betay*cos(conj(theta))*sin(conj(fi))*sin(fi)*sin(theta))/m + (betay*sin(conj(fi))*sin(conj(theta))*cos(theta)*sin(fi))/m) + vx*((betay*sin(fi)*sin(theta)*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta))))/m - (betaz*cos(fi)*sin(theta)*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta))))/m - (betax*cos(conj(psi))*cos(conj(theta))*cos(theta))/m + (betax*cos(conj(psi))*sin(conj(theta))*sin(theta))/m + (betaz*cos(conj(fi))*cos(conj(psi))*cos(conj(theta))*cos(fi)*cos(theta))/m + (betay*cos(conj(psi))*cos(conj(theta))*sin(conj(fi))*cos(theta)*sin(fi))/m))  , abs(vz)*(vy*((betay*cos(theta)*sin(fi)*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta))))/m - (betaz*cos(fi)*cos(theta)*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta))))/m + (betax*cos(conj(psi))*cos(conj(theta))*sin(theta))/m) - vx*((betaz*cos(fi)*cos(theta)*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m - (betay*cos(theta)*sin(fi)*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m + (betax*cos(conj(theta))*sin(conj(psi))*sin(theta))/m))                  ] ;
%v/lbd
   
A11_11 = sign(vx)*(vy*((betay*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m + (betaz*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m - (betax*cos(conj(theta))*sin(conj(psi))*cos(psi)*cos(theta))/m) - vx*((betaz*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta)))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m + (betay*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta)))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m + (betax*cos(conj(psi))*cos(conj(theta))*cos(psi)*cos(theta))/m) + vz*((betay*cos(conj(theta))*sin(conj(fi))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m - (betaz*cos(conj(fi))*cos(conj(theta))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m + (betax*sin(conj(theta))*cos(psi)*cos(theta))/m)) - abs(vx)*((betaz*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta)))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m + (betay*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta)))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m + (betax*cos(conj(psi))*cos(conj(theta))*cos(psi)*cos(theta))/m)   ;
A11_12 = wz + abs(vx)*((betay*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m + (betaz*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m - (betax*cos(conj(theta))*sin(conj(psi))*cos(psi)*cos(theta))/m);
A11_13 = abs(vx)*((betay*cos(conj(theta))*sin(conj(fi))*(cos(fi)*sin(psi) - cos(psi)*sin(fi)*sin(theta)))/m - (betaz*cos(conj(fi))*cos(conj(theta))*(sin(fi)*sin(psi) + cos(fi)*cos(psi)*sin(theta)))/m + (betax*sin(conj(theta))*cos(psi)*cos(theta))/m) - wy; 
A11_21 = abs(vy)*((betay*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta)))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m + (betaz*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta)))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m - (betax*cos(conj(psi))*cos(conj(theta))*cos(theta)*sin(psi))/m) - wz  ; 
A11_22 = sign(vy)*(vx*((betay*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta)))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m + (betaz*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta)))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m - (betax*cos(conj(psi))*cos(conj(theta))*cos(theta)*sin(psi))/m) - vy*((betay*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m + (betaz*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m + (betax*cos(conj(theta))*sin(conj(psi))*cos(theta)*sin(psi))/m) + vz*((betaz*cos(conj(fi))*cos(conj(theta))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m - (betay*cos(conj(theta))*sin(conj(fi))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m + (betax*sin(conj(theta))*cos(theta)*sin(psi))/m)) - abs(vy)*((betay*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m + (betaz*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta)))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m + (betax*cos(conj(theta))*sin(conj(psi))*cos(theta)*sin(psi))/m); 
A11_23 = wx + abs(vy)*((betaz*cos(conj(fi))*cos(conj(theta))*(cos(psi)*sin(fi) - cos(fi)*sin(psi)*sin(theta)))/m - (betay*cos(conj(theta))*sin(conj(fi))*(cos(fi)*cos(psi) + sin(fi)*sin(psi)*sin(theta)))/m + (betax*sin(conj(theta))*cos(theta)*sin(psi))/m); 
A11_31 = wy + abs(vz)*((betay*cos(theta)*sin(fi)*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta))))/m - (betaz*cos(fi)*cos(theta)*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta))))/m + (betax*cos(conj(psi))*cos(conj(theta))*sin(theta))/m); 
A11_32 = abs(vz)*((betaz*cos(fi)*cos(theta)*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m - (betay*cos(theta)*sin(fi)*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m + (betax*cos(conj(theta))*sin(conj(psi))*sin(theta))/m) - wx; 
A11_33 = sign(vz)*(vx*((betay*cos(theta)*sin(fi)*(cos(conj(fi))*sin(conj(psi)) - cos(conj(psi))*sin(conj(fi))*sin(conj(theta))))/m - (betaz*cos(fi)*cos(theta)*(sin(conj(fi))*sin(conj(psi)) + cos(conj(fi))*cos(conj(psi))*sin(conj(theta))))/m + (betax*cos(conj(psi))*cos(conj(theta))*sin(theta))/m) + vy*((betaz*cos(fi)*cos(theta)*(cos(conj(psi))*sin(conj(fi)) - cos(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m - (betay*cos(theta)*sin(fi)*(cos(conj(fi))*cos(conj(psi)) + sin(conj(fi))*sin(conj(psi))*sin(conj(theta))))/m + (betax*cos(conj(theta))*sin(conj(psi))*sin(theta))/m) - vz*((betax*sin(conj(theta))*sin(theta))/m + (betaz*cos(conj(fi))*cos(conj(theta))*cos(fi)*cos(theta))/m + (betay*cos(conj(theta))*sin(conj(fi))*cos(theta)*sin(fi))/m)) - abs(vz)*((betax*sin(conj(theta))*sin(theta))/m + (betaz*cos(conj(fi))*cos(conj(theta))*cos(fi)*cos(theta))/m + (betay*cos(conj(theta))*sin(conj(fi))*cos(theta)*sin(fi))/m); 

A11 =   [A11_11 , A11_12 , A11_13 ;
         A11_21 , A11_22 , A11_23 ;
         A11_31 , A11_32 , A11_13 ]; 
%v/v
 
A12 = [  0 , -vz  ,   vy  ;
        vz ,   0  ,  -vx ;
       -vy ,  vx  ,    0 ] ;
%v/om

A13 = zeros (3,3) ;
%om/p
 
A14 = zeros (3,3) ;
%om/lbd
 
A15 = zeros (3,3) ;
%om/v

A16 = [         0              ,  (Jy*wz)/Jx - (Jz*wz)/Jx   , (Jy*wy)/Jx - (Jz*wy)/Jx  ;
       (Jz*wz)/Jy - (Jx*wz)/Jy ,   0                        , (Jz*wx)/Jy - (Jx*wx)/Jy  ;
       (Jx*wy)/Jz - (Jy*wy)/Jz ,  (Jx*wx)/Jz - (Jy*wx)/Jz   ,             0            ] ;
%om/om
  
A = [A1  , A3  , A2  , A4  ;
     A9  , A11 , A10 , A12 ;
     A5  , A7  , A6  , A8  ;
     A13 , A15 , A14 , A16 ] 
  
B1 = zeros (3,1) ; %p/T

B2 = zeros (3,1) ; %p/nx

B3 = zeros (3,1) ; %p/ny

B4 = zeros (3,1) ; %p/nz

B5 = 1/m*[0 ; 0 ; 1] ; %v/T

B6 = zeros (3,1) ; %v/nx

B7 = zeros (3,1) ; %v/ny

B8 = zeros (3,1) ; %v/nz

B9 =  zeros (3,1) ; %lbd/T

B10 =  zeros (3,1) ; %lbd/nx

B11 =  zeros (3,1) ; %lbd/ny

B12 =  zeros (3,1) ; %lbd/nz

B13 = zeros(3,1) ; %om/T

B14 = [1/Jx ; 0 ; 0 ] ;  %om/nx

B15 = [0  ; 1/Jy ; 0 ] ; %om/ny
 
B16 = [ 0 ; 0 ; 1/Jz ] ;  %om/nz

B = [B1  , B2  , B3  , B4  ;  
     B5  , B6  , B7  , B8  ;
     B9  , B10 , B11 , B12 ;
     B13 , B14 , B15 , B16 ] 

C = eye(12);
 
D = zeros(12,4);
 
sys = ss(A,B,C,D);
 
u_L = [ T1 ; T2 ; T3 ; T4 ]*(t>=0);
 
%simulate linear system:
y_L = lsim(sys,u_L,t,xop)';

% Calcule os autovalores e autovetores
[autovetores, autovalores] = eig(A);

% Exiba os autovalores
disp('Autovalores:');
disp(diag(autovalores));

% Exiba os autovetores
disp('Autovetores:');
disp(autovetores);

% Analise a estabilidade do sistema
partes_reais = real(diag(autovalores));
if all(partes_reais < 0)
    disp('O sistema é estável.');
elseif any(partes_reais > 0)
    disp('O sistema é instável.');
else
    disp('O sistema é marginalmente estável.');
end

% test controlability, observability and stability
[V,D1,W] = eig(A)
[Vj,Jor] = jordan(A);
Jor

mode_ctrl = W'*B
n1_unstable_modes = rank(ctrb(A,B))-12;
if n1_unstable_modes > 0, disp('Linearized system is not controlable.'); else disp('Linearized system is controlable.'); end

mode_obs = C*V
n1_unobservable_modes = rank(obsv(A,C))-12;
if n1_unobservable_modes > 0, disp('Linearized system is not observable.'); else disp('Linearized system is observable.'); end

Gss = ss(A,B,C,D);
G = zpk(tf(Gss))

%py,fi
G0 = G(2,2);
G1 = G(7,2)^-1;
G2 = G(2,2)*G(7,2)^-1

%px,theta
a =  G(1,3);
b = G(8,3);
G3 = G(1,3)*G(8,3)^-1


figure(2);
set(groot,'defaultLineLineWidth',0.5)
rlocus(G(3,1));
title('Root Locus G_{T,pz}');
grid on


figure(3);
set(groot,'defaultLineLineWidth',0.5)
rlocus(G3);
title('Root Locus G_{\theta,px}');
grid on