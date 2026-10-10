clc;
clear all; %#ok<*CLALL>
close all;
global m1 l1 m2 l2 m3 l3 m4 l4 k1 k2 k3
m1 = 3.5; m2 = 3; m3 = 2.5; m4 = 2;
l1 = 2; l2 = 1.5; l3 = 1; l4 = 0.5;
k1 = 3; k2 = 4; k3 = 2;

%_______________Choosing the method of numerical solution_________________%
choice1 = menu('Please Choose the Method','Integrated Multipliers','Augmented','Elimination'...
    ,'Embedded','Modified Lagrange','Greenwood','Gibbs-Appel','Kane');
IM = 0; AU = 0; EL = 0; EM = 0; ML = 0; GR = 0; Kane = 0; Gibss = 0; Anim=0;
if choice1 == 1
    IM = 1;   %1 to use Integrated Multipliers method
end
if choice1 == 2
    AU = 1;   %1 to use Augmented method
end
if choice1 == 3
    EL = 1;   %1 to use Elimination method
end
if choice1 == 4
    EM = 1;   %1 to use Embedding method
end
if choice1 == 5
    ML = 1;   %1 to use Modified Lagrange method
end
if choice1 == 6
   GR = 1;   %1 to use Greenwood method
end
if choice1 == 7
   Gibss = 1; 
end
if choice1 == 8
   Kane = 1;
end
choice2 = menu('Please choose whether to show animation or plots','Animation','Plots');

if choice2 == 1
    Anim = 1; %1 to show animation instead of plots
end
if choice2 == 2
    Anim = 0; %1 to show animation instead of plots
end
simulation_time = [0:0.01:5]; %#ok<*NBRAK>

%__________________Integrated Multipliers Method__________________________%
if IM==1
z0 = [1;1.5; pi/12;pi/4;pi/10;pi/6;0;0;0;0;0;0;0;0;0;0];
options=odeset('maxstep',0.001);
tic;
[t,z] = ode45(@IntMulDyn,simulation_time,z0, options);
toc
q1_IM = z(:,1); q2_IM = z(:,2); q3_IM = z(:,3); q4_IM = z(:,4); q5_IM = z(:,5); q6_IM = z(:,6);
dq1_IM = z(:,7); dq2_IM = z(:,8); dq3_IM = z(:,9); dq4_IM = z(:,10); dq5_IM = z(:,11); dq6_IM = z(:,12);
lamda1_IM = z(:,13); lamda2_IM = z(:,14); lamda3_IM = z(:,15); lamda4_IM = z(:,16);

% Plotting Generalized Coordinates
if Anim == 0
figure;
hold on;
plot(t,q1_IM ,'r-','linewidth',4)
plot(t,q2_IM ,'b-','linewidth',4)
legend('x(m)','y(m)')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');
figure;
hold on;
plot(t,q3_IM,'b-','linewidth',4)
plot(t,q4_IM,'g-','linewidth',4)
plot(t,q5_IM,'k-','linewidth',4)
plot(t,q6_IM,'r-','linewidth',4)
legend('\theta_1','\theta_2','\theta_3','\theta_4')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');

% Plotting Constraints
Cons_IM = constraints(dq1_IM,dq2_IM,dq3_IM,dq4_IM,dq5_IM,dq6_IM,l1,l2,l3,l4,q3_IM,q4_IM,q5_IM,q6_IM);
figure;
hold on;
plot(t,Cons_IM(1:length(simulation_time)) ,'r-','linewidth',4)
plot(t,Cons_IM(1+length(simulation_time):2*length(simulation_time)) ,'b-','linewidth',4)
plot(t,Cons_IM(1+2*length(simulation_time):3*length(simulation_time)) ,'g-','linewidth',4)
plot(t,Cons_IM(1+3*length(simulation_time):4*length(simulation_time)) ,'k-','linewidth',4)
legend('f1','f2','f3','f4')
set(gca,'fontsize',14,'fontweight','bold');
xlabel('Time(s)');
ylabel('Error of Constraints','fontsize',14,'fontweight','bold');

% Plotting Energy Error
E_IM = Energy(dq1_IM,dq2_IM,dq3_IM,dq4_IM,dq5_IM,dq6_IM,k1,k2,k3,l1,l2,l3,l4,m1,m2,m3,m4,q3_IM,q4_IM,q5_IM,q6_IM);
figure;
plot(t,100*(E_IM-E_IM(1))/E_IM(1),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Relative Error of  Mechanical Energy % ','fontsize',14,'fontweight','bold');
end
end


%__________________________Augmented Method_______________________________%
if AU==1
z0 = [1;1.5; pi/12;pi/4;pi/10;pi/6;0;0;0;0;0;0;0;0;0;0];
options=odeset('maxstep',0.001);
tic;
[t,z] = ode45(@AugmntDyn,simulation_time,z0, options);
toc
q1_Au = z(:,1); q2_Au = z(:,2); q3_Au = z(:,3); q4_Au = z(:,4); q5_Au = z(:,5); q6_Au = z(:,6);
dq1_Au = z(:,7); dq2_Au = z(:,8); dq3_Au = z(:,9); dq4_Au = z(:,10); dq5_Au = z(:,11); dq6_Au = z(:,12);
lamda1_Au = z(:,13); lamda2_Au = z(:,14); lamda3_Au = z(:,15); lamda4_Au = z(:,16);

if Anim == 0
% Plotting Generalized Coordinates
figure;
hold on;
plot(t,q1_Au ,'r-','linewidth',4)
plot(t,q2_Au ,'b-','linewidth',4)
legend('x(m)','y(m)')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');
figure;
hold on;
plot(t,q3_Au,'b-','linewidth',4)
plot(t,q4_Au,'g-','linewidth',4)
plot(t,q5_Au,'k-','linewidth',4)
plot(t,q6_Au,'r-','linewidth',4)
legend('\theta_1','\theta_2','\theta_3','\theta_4')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');

% Plotting Constraints
Cons_Au = constraints(dq1_Au,dq2_Au,dq3_Au,dq4_Au,dq5_Au,dq6_Au,l1,l2,l3,l4,q3_Au,q4_Au,q5_Au,q6_Au);
figure;
hold on;
plot(t,Cons_Au(1:length(simulation_time)) ,'r-','linewidth',4)
plot(t,Cons_Au(1+length(simulation_time):2*length(simulation_time)) ,'b-','linewidth',4)
plot(t,Cons_Au(1+2*length(simulation_time):3*length(simulation_time)) ,'g-','linewidth',4)
plot(t,Cons_Au(1+3*length(simulation_time):4*length(simulation_time)) ,'k-','linewidth',4)
legend('f1','f2','f3','f4')
set(gca,'fontsize',14,'fontweight','bold');
xlabel('Time(s)');
ylabel('Error of Constraints','fontsize',14,'fontweight','bold');

% Plotting Energy Error
E_Au = Energy(dq1_Au,dq2_Au,dq3_Au,dq4_Au,dq5_Au,dq6_Au,k1,k2,k3,l1,l2,l3,l4,m1,m2,m3,m4,q3_Au,q4_Au,q5_Au,q6_Au);
figure;
plot(t,100*(E_Au-E_Au(1))/E_Au(1),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Relative Error of  Mechanical Energy % ','fontsize',14,'fontweight','bold');
end
end

%___________________________Elimination Method____________________________%
if EL==1
z0 = [1;1.5; pi/12;pi/4;pi/10;pi/6;0;0;0;0;0;0];
options=odeset('maxstep',0.001);
tic;
[t,z] = ode45(@ElimDyn,simulation_time,z0, options);
toc
q1_El = z(:,1); q2_El = z(:,2); q3_El = z(:,3); q4_El = z(:,4); q5_El = z(:,5); q6_El = z(:,6);
dq1_El = z(:,7); dq2_El = z(:,8); dq3_El = z(:,9); dq4_El = z(:,10); dq5_El = z(:,11); dq6_El = z(:,12);

% Plotting Generalized Coordinates
if Anim == 0
figure;
hold on;
plot(t,q1_El ,'r-','linewidth',4)
plot(t,q2_El ,'b-','linewidth',4)
legend('x(m)','y(m)')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');
figure;
hold on;
plot(t,q3_El,'b-','linewidth',4)
plot(t,q4_El,'g-','linewidth',4)
plot(t,q5_El,'k-','linewidth',4)
plot(t,q6_El,'r-','linewidth',4)
legend('\theta_1','\theta_2','\theta_3','\theta_4')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');

% Plotting Constraints
Cons_El = constraints(dq1_El,dq2_El,dq3_El,dq4_El,dq5_El,dq6_El,l1,l2,l3,l4,q3_El,q4_El,q5_El,q6_El);
figure;
hold on;
plot(t,Cons_El(1:length(simulation_time)) ,'r-','linewidth',4)
plot(t,Cons_El(1+length(simulation_time):2*length(simulation_time)) ,'b-','linewidth',4)
plot(t,Cons_El(1+2*length(simulation_time):3*length(simulation_time)) ,'g-','linewidth',4)
plot(t,Cons_El(1+3*length(simulation_time):4*length(simulation_time)) ,'k-','linewidth',4)
legend('f1','f2','f3','f4')
set(gca,'fontsize',14,'fontweight','bold');
xlabel('Time(s)');
ylabel('Error of Constraints','fontsize',14,'fontweight','bold');

% Plotting Energy Error
E_El = Energy(dq1_El,dq2_El,dq3_El,dq4_El,dq5_El,dq6_El,k1,k2,k3,l1,l2,l3,l4,m1,m2,m3,m4,q3_El,q4_El,q5_El,q6_El);
figure;
plot(t,100*(E_El-E_El(1))/E_El(1),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Relative Error of  Mechanical Energy % ','fontsize',14,'fontweight','bold');
end
end


%____________________________Greenwood Method_____________________________%
if GR==1
z0 = [1;1.5; pi/12;pi/4;pi/10;pi/6;0;0;0;0;0;0];
options=odeset('maxstep',0.001);
tic;
[t,z] = ode45(@GreenDyn,simulation_time,z0, options);
toc
q1_Gr = z(:,1); q2_Gr = z(:,2); q3_Gr = z(:,3); q4_Gr = z(:,4); q5_Gr = z(:,5); q6_Gr = z(:,6);
dq1_Gr = z(:,7); dq2_Gr = z(:,8); dq3_Gr = z(:,9); dq4_Gr = z(:,10); dq5_Gr = z(:,11); dq6_Gr = z(:,12);

% Plotting Generalized Coordinates
if Anim == 0
figure;
hold on;
plot(t,q1_Gr ,'r-','linewidth',4)
plot(t,q2_Gr ,'b-','linewidth',4)
legend('x(m)','y(m)')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');
figure;
hold on;
plot(t,q3_Gr,'b-','linewidth',4)
plot(t,q4_Gr,'g-','linewidth',4)
plot(t,q5_Gr,'k-','linewidth',4)
plot(t,q6_Gr,'r-','linewidth',4)
legend('\theta_1','\theta_2','\theta_3','\theta_4')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');

% Plotting Constraints
Cons_Gr = constraints(dq1_Gr,dq2_Gr,dq3_Gr,dq4_Gr,dq5_Gr,dq6_Gr,l1,l2,l3,l4,q3_Gr,q4_Gr,q5_Gr,q6_Gr);
figure;
hold on;
plot(t,Cons_Gr(1:length(simulation_time)) ,'r-','linewidth',4)
plot(t,Cons_Gr(1+length(simulation_time):2*length(simulation_time)) ,'b-','linewidth',4)
plot(t,Cons_Gr(1+2*length(simulation_time):3*length(simulation_time)) ,'g-','linewidth',4)
plot(t,Cons_Gr(1+3*length(simulation_time):4*length(simulation_time)) ,'k-','linewidth',4)
legend('f1','f2','f3','f4')
set(gca,'fontsize',14,'fontweight','bold');
xlabel('Time(s)');
ylabel('Error of Constraints','fontsize',14,'fontweight','bold');

% Plotting Energy Error
E_Gr = Energy(dq1_Gr,dq2_Gr,dq3_Gr,dq4_Gr,dq5_Gr,dq6_Gr,k1,k2,k3,l1,l2,l3,l4,m1,m2,m3,m4,q3_Gr,q4_Gr,q5_Gr,q6_Gr);
figure;
plot(t,100*(E_Gr-E_Gr(1))/E_Gr(1),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Relative Error of  Mechanical Energy % ','fontsize',14,'fontweight','bold');
end
end

%________________________Modified Lagrange Method_________________________%
if ML==1
z0 = [1;1.5; pi/12;pi/4;pi/10;pi/6;0;0];
options=odeset('maxstep',0.001);
tic
[t,z] = ode45(@ModlagDyn,simulation_time,z0, options);
toc;
q1_Ml = z(:,1); q2_Ml = z(:,2); q3_Ml = z(:,3); q4_Ml = z(:,4); q5_Ml = z(:,5); q6_Ml = z(:,6);
u1 = z(:,7); u2 = z(:,8);
qdot_Ml = zeros(6,length(simulation_time));
for i=1:length(simulation_time)
u = [u1(i);u2(i)];
W = W_Ml(l1,l2,l3,l4,q3_Ml(i),q4_Ml(i),q5_Ml(i),q6_Ml(i));
X = X_Ml;
qdot_Ml(:,i) = W*u + X;
end

% Plotting Generalized Coordinates
if Anim == 0
figure;
hold on;
plot(t,q1_Ml ,'r-','linewidth',4)
plot(t,q2_Ml ,'b-','linewidth',4)
legend('x(m)','y(m)')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');
figure;
hold on;
plot(t,q3_Ml,'b-','linewidth',4)
plot(t,q4_Ml,'g-','linewidth',4)
plot(t,q5_Ml,'k-','linewidth',4)
plot(t,q6_Ml,'r-','linewidth',4)
legend('\theta_1','\theta_2','\theta_3','\theta_4')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');

% Plotting Constraints
Cons_Ml = constraints(qdot_Ml(1,:)',qdot_Ml(2,:)',qdot_Ml(3,:)',qdot_Ml(4,:)',qdot_Ml(5,:)',qdot_Ml(6,:)',l1,l2,l3,l4,q3_Ml,q4_Ml,q5_Ml,q6_Ml);
figure;
hold on;
plot(t,Cons_Ml(1:length(simulation_time)) ,'r-','linewidth',4)
plot(t,Cons_Ml(1+length(simulation_time):2*length(simulation_time)) ,'b-','linewidth',4)
plot(t,Cons_Ml(1+2*length(simulation_time):3*length(simulation_time)) ,'g-','linewidth',4)
plot(t,Cons_Ml(1+3*length(simulation_time):4*length(simulation_time)) ,'k-','linewidth',4)
legend('f1','f2','f3','f4')
set(gca,'fontsize',14,'fontweight','bold');
xlabel('Time(s)');
ylabel('Error of Constraints','fontsize',14,'fontweight','bold');

% Plotting Energy Error
E_Ml = Energy(qdot_Ml(1,:)',qdot_Ml(2,:)',qdot_Ml(3,:)',qdot_Ml(4,:)',qdot_Ml(5,:)',qdot_Ml(6,:)',k1,k2,k3,l1,l2,l3,l4,m1,m2,m3,m4,q3_Ml,q4_Ml,q5_Ml,q6_Ml);
figure;
plot(t,100*(E_Ml-E_Ml(1))/E_Ml(1),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Relative Error of  Mechanical Energy % ','fontsize',14,'fontweight','bold');
end
end

%____________________________Embedding Method_____________________________%
if EM==1
z0 = [1;1.5; pi/12;pi/4;pi/10;pi/6;0;0];
options=odeset('maxstep',0.001);
tic;
[t,z] = ode45(@EmbdDyn,simulation_time,z0, options);
toc
q1_Em = z(:,1); q2_Em = z(:,2); q3_Em = z(:,3); q4_Em = z(:,4); q5_Em = z(:,5); q6_Em = z(:,6);
dq1 = z(:,7); dq6 = z(:,8);
qdot_Em = zeros(6,length(simulation_time));
for i=1:length(simulation_time)
dqI = [dq1(i);dq6(i)];
Wbar = Wbar_Embd(l1,l2,l3,l4,q3_Em(i),q4_Em(i),q5_Em(i),q6_Em(i));
Xbar = Xbar_Embd;
qdot_Em(:,i) = Wbar*dqI + Xbar;
end

% Plotting Generalized Coordinates
if Anim == 0
figure;
hold on;
plot(t,q1_Em ,'r-','linewidth',4)
plot(t,q2_Em ,'b-','linewidth',4)
legend('x(m)','y(m)')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');
figure;
hold on;
plot(t,q3_Em,'b-','linewidth',4)
plot(t,q4_Em,'g-','linewidth',4)
plot(t,q5_Em,'k-','linewidth',4)
plot(t,q6_Em,'r-','linewidth',4)
legend('\theta_1','\theta_2','\theta_3','\theta_4')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');

% Plotting Constraints
Cons_Em = constraints(qdot_Em(1,:)',qdot_Em(2,:)',qdot_Em(3,:)',qdot_Em(4,:)',qdot_Em(5,:)',qdot_Em(6,:)',l1,l2,l3,l4,q3_Em,q4_Em,q5_Em,q6_Em);
figure;
hold on;
plot(t,Cons_Em(1:length(simulation_time)) ,'r-','linewidth',4)
plot(t,Cons_Em(1+length(simulation_time):2*length(simulation_time)) ,'b-','linewidth',4)
plot(t,Cons_Em(1+2*length(simulation_time):3*length(simulation_time)) ,'g-','linewidth',4)
plot(t,Cons_Em(1+3*length(simulation_time):4*length(simulation_time)) ,'k-','linewidth',4)
legend('f1','f2','f3','f4')
set(gca,'fontsize',14,'fontweight','bold');
xlabel('Time(s)');
ylabel('Error of Constraints','fontsize',14,'fontweight','bold');

% Plotting Energy Error
E_Em = Energy(qdot_Em(1,:)',qdot_Em(2,:)',qdot_Em(3,:)',qdot_Em(4,:)',qdot_Em(5,:)',qdot_Em(6,:)',k1,k2,k3,l1,l2,l3,l4,m1,m2,m3,m4,q3_Em,q4_Em,q5_Em,q6_Em);
figure;
plot(t,100*(E_Em-E_Em(1))/E_Em(1),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Relative Error of  Mechanical Energy % ','fontsize',14,'fontweight','bold');
end
end





%___________________________Gibss-Appel Method____________________________%
if Gibss == 1
z0 = [1;1.5; pi/12;pi/4;pi/10;pi/6;0;0];
options=odeset('maxstep',0.001);
tic;
[t,z] = ode45(@GibssDyn,simulation_time,z0, options);
toc
q1_gb = z(:,1); q2_gb = z(:,2); q3_gb = z(:,3); q4_gb = z(:,4); q5_gb = z(:,5); q6_gb = z(:,6);
u1 = z(:,7); u2 = z(:,8);
qdot_gb = zeros(6,length(simulation_time));
for i=1:length(simulation_time)
u = [u1(i);u2(i)];
W = W_Ml(l1,l2,l3,l4,q3_gb(i),q4_gb(i),q5_gb(i),q6_gb(i));
X = X_Ml;
qdot_gb(:,i) = W*u + X;
end

% Plotting Generalized Coordinates
if Anim == 0
figure;
hold on;
plot(t,q1_gb ,'r-','linewidth',4)
plot(t,q2_gb ,'b-','linewidth',4)
legend('x(m)','y(m)')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');
figure;
hold on;
plot(t,q3_gb,'b-','linewidth',4)
plot(t,q4_gb,'g-','linewidth',4)
plot(t,q5_gb,'k-','linewidth',4)
plot(t,q6_gb,'r-','linewidth',4)
legend('\theta_1','\theta_2','\theta_3','\theta_4')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');

% Plotting Constraints
Cons_gb = constraints(qdot_gb(1,:)',qdot_gb(2,:)',qdot_gb(3,:)',qdot_gb(4,:)',qdot_gb(5,:)',qdot_gb(6,:)',l1,l2,l3,l4,q3_gb,q4_gb,q5_gb,q6_gb);
figure;
hold on;
plot(t,Cons_gb(1:length(simulation_time)) ,'r-','linewidth',4)
plot(t,Cons_gb(1+length(simulation_time):2*length(simulation_time)) ,'b-','linewidth',4)
plot(t,Cons_gb(1+2*length(simulation_time):3*length(simulation_time)) ,'g-','linewidth',4)
plot(t,Cons_gb(1+3*length(simulation_time):4*length(simulation_time)) ,'k-','linewidth',4)
legend('f1','f2','f3','f4')
set(gca,'fontsize',14,'fontweight','bold');
xlabel('Time(s)');
ylabel('Error of Constraints','fontsize',14,'fontweight','bold');

% Plotting Energy Error
E_gb = Energy(qdot_gb(1,:)',qdot_gb(2,:)',qdot_gb(3,:)',qdot_gb(4,:)',qdot_gb(5,:)',qdot_gb(6,:)',k1,k2,k3,l1,l2,l3,l4,m1,m2,m3,m4,q3_gb,q4_gb,q5_gb,q6_gb);
figure;
plot(t,100*(E_gb-E_gb(1))/E_gb(1),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Relative Error of  Mechanical Energy % ','fontsize',14,'fontweight','bold');
end
end





%___________________________Kane Method____________________________%
if Kane == 1
z0 = [1;1.5; pi/12;pi/4;pi/10;pi/6;0;0];
options=odeset('maxstep',0.001);
tic;
[t,z] = ode45(@KaneDyn,simulation_time,z0, options);
toc
q1_kn = z(:,1); q2_kn = z(:,2); q3_kn = z(:,3); q4_kn = z(:,4); q5_kn = z(:,5); q6_kn = z(:,6);
u1 = z(:,7); u2 = z(:,8);
qdot_kn = zeros(6,length(simulation_time));
for i=1:length(simulation_time)
u = [u1(i);u2(i)];
W = W_Ml(l1,l2,l3,l4,q3_kn(i),q4_kn(i),q5_kn(i),q6_kn(i));
X = X_Ml;
qdot_kn(:,i) = W*u + X;
end

% Plotting Generalized Coordinates
if Anim == 0
figure;
hold on;
plot(t,q1_kn ,'r-','linewidth',4)
plot(t,q2_kn ,'b-','linewidth',4)
legend('x(m)','y(m)')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');
figure;
hold on;
plot(t,q3_kn,'b-','linewidth',4)
plot(t,q4_kn,'g-','linewidth',4)
plot(t,q5_kn,'k-','linewidth',4)
plot(t,q6_kn,'r-','linewidth',4)
legend('\theta_1','\theta_2','\theta_3','\theta_4')
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');


% Plotting Constraints
Cons_kn = constraints(qdot_kn(1,:)',qdot_kn(2,:)',qdot_kn(3,:)',qdot_kn(4,:)',qdot_kn(5,:)',qdot_kn(6,:)',l1,l2,l3,l4,q3_kn,q4_kn,q5_kn,q6_kn);
figure;
hold on;
plot(t,Cons_kn(1:length(simulation_time)) ,'r-','linewidth',4)
plot(t,Cons_kn(1+length(simulation_time):2*length(simulation_time)) ,'b-','linewidth',4)
plot(t,Cons_kn(1+2*length(simulation_time):3*length(simulation_time)) ,'g-','linewidth',4)
plot(t,Cons_kn(1+3*length(simulation_time):4*length(simulation_time)) ,'k-','linewidth',4)
legend('f1','f2','f3','f4')
set(gca,'fontsize',14,'fontweight','bold');
xlabel('Time(s)');
ylabel('Error of Constraints','fontsize',14,'fontweight','bold');

% Plotting Energy Error
E_kn = Energy(qdot_kn(1,:)',qdot_kn(2,:)',qdot_kn(3,:)',qdot_kn(4,:)',qdot_kn(5,:)',qdot_kn(6,:)',k1,k2,k3,l1,l2,l3,l4,m1,m2,m3,m4,q3_kn,q4_kn,q5_kn,q6_kn);
figure;
plot(t,100*(E_kn-E_kn(1))/E_kn(1),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Relative Error of  Mechanical Energy % ','fontsize',14,'fontweight','bold');
end
end







%_______________________________Animation_________________________________%
if Anim ==1
    if IM == 1
q1 = q1_IM; q2 = q2_IM; q3 = q3_IM; q4 = q4_IM; q5 = q5_IM; q6 = q6_IM;  
    end
    if AU == 1
q1 = q1_Au; q2 = q2_Au; q3 = q3_Au; q4 = q4_Au; q5 = q5_Au; q6 = q6_Au;  
    end
    if EL == 1
q1 = q1_El; q2 = q2_El; q3 = q3_El; q4 = q4_El; q5 = q5_El; q6 = q6_El;  
    end
    if EM == 1
q1 = q1_Em; q2 = q2_Em; q3 = q3_Em; q4 = q4_Em; q5 = q5_Em; q6 = q6_Em;  
    end
    if ML == 1
q1 = q1_Ml; q2 = q2_Ml; q3 = q3_Ml; q4 = q4_Ml; q5 = q5_Ml; q6 = q6_Ml;  
    end
    if GR == 1
q1 = q1_Gr; q2 = q2_Gr; q3 = q3_Gr; q4 = q4_Gr; q5 = q5_Gr; q6 = q6_Gr;  
    end
    if Gibss == 1
q1 = q1_gb; q2 = q2_gb; q3 = q3_gb; q4 = q4_gb; q5 = q5_gb; q6 = q6_gb;  
    end
    if Kane == 1
q1 = q1_kn; q2 = q2_kn; q3 = q3_kn; q4 = q4_kn; q5 = q5_kn; q6 = q6_kn;  
    end
figure;
for i=1:5:length(simulation_time)
    xG1 = q1(i); yG1 = q2(i);
    xS = xG1 - (l1/2)*cos(q3(i)); yS = yG1 - (l1/2)*sin(q3(i));     % S is the start point of the first link
    xA = xG1 + (l1/2)*cos(q3(i)); yA = yG1 + (l1/2)*sin(q3(i)); 
    xB = xA + l2*cos(q4(i)); yB = yA + l2*sin(q4(i));
    xC = xB + l3*cos(q5(i)); yC = yB + l3*sin(q5(i));
    xD = xC + l4*cos(q6(i)); yD = yC + l4*sin(q6(i));     % D is the end of the last link
    
    
    v1x = [xS; xA]; v1y = [yS; yA];
    v2x = [xA; xB]; v2y = [yA; yB];
    v3x = [xB; xC]; v3y = [yB; yC];
    v4x = [xC; xD]; v4y = [yC; yD];
    
    
    plot(v1x,v1y,'k-','linewidth',6)
    hold on
    plot(v2x,v2y,'k-','linewidth',6)
    plot(v3x,v3y,'k-','linewidth',6)
    plot(v4x,v4y,'k-','linewidth',6)
    
    plot(xS,yS,'bO','linewidth',4)
    plot(xA,yA,'rO','linewidth',4)
    plot(xB,yB,'rO','linewidth',4)
    plot(xC,yC,'rO','linewidth',4)
    plot(xD,yD,'bO','linewidth',4)
    grid on
    axis equal
    axis([-l2  (l4+l3+l2+l1) -l2 (l4+l3+l2+l1)])
    
    hold off
    pause(0.001)
end
end