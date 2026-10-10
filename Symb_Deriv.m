clc;
clear all; %#ok<*CLALL>
close all;

syms q1 q2 q3 q4 q5 q6 real     %x, y, theta1, theta2, theta3, theta4
syms dq1 dq2 dq3 dq4 dq5 dq6 real
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 real
syms U1 U2 dU1 dU2 real
syms t
syms m1 l1 m2 l2 m3 l3 m4 l4 k1 k2 k3 real

q = [q1;q2;q3;q4;q5;q6];
dq = [dq1;dq2;dq3;dq4;dq5;dq6];
ddq = [ddq1;ddq2;ddq3;ddq4;ddq5;ddq6];

% Angular Velocities
w1 = [0;0;dq3];
w2 = [0;0;dq4];
w3 = [0;0;dq5];
w4 = [0;0;dq6];

% Position Vectors
r_A_G1 = (l1/2)*[cos(q3);sin(q3);0];
r_G2_A = (l2/2)*[cos(q4);sin(q4);0];
r_B_G2 = (l2/2)*[cos(q4);sin(q4);0];
r_G3_B = (l3/2)*[cos(q5);sin(q5);0];
r_C_G3 = (l3/2)*[cos(q5);sin(q5);0];
r_G4_C = (l4/2)*[cos(q6);sin(q6);0];

% Transitional Velocities
v_G1 = [dq1;dq2;0];
v_A = v_G1 + cross(w1,r_A_G1);
v_G2 = v_A + cross(w2,r_G2_A);
v_B = v_G2 + cross(w2,r_B_G2);
v_G3 = v_B + cross(w3,r_G3_B);
v_C = v_G3 + cross(w3,r_C_G3);
v_G4 = v_C + cross(w4,r_G4_C);

% Constraints' Directions
e1 = [-sin(q3);cos(q3);0];     %e1 is normal to the link1's direction
e2 = [-sin(q4);cos(q4);0];
e3 = [-sin(q5);cos(q5);0];
e4 = [-sin(q6);cos(q6);0];

% Constraints
f1 = dot(v_G1,e1);
f2 = dot(v_G2,e2);
f3 = dot(v_G3,e3);
f4 = dot(v_G4,e4);
f = [f1;f2;f3;f4];
% matlabFunction(f, 'File', 'constraints');

a = jacobian(f,dq);
a = simplify(a);
b = f - a*dq;
b = simplify(b);

% Derivatives of a and b
adot = zeros(size(a));
bdot = zeros(size(b));
for i=1:length(q)
    adot = adot + diff(a,q(i))*dq(i);
    bdot = bdot + diff(b,q(i))*dq(i);
end
adot = simplify(adot);
bdot = simplify(bdot);

% Kinetic Energy
T1 = 0.5*m1*(v_G1.'*v_G1) + 0.5*(m1*l1*l1/12)*(w1.'*w1);
T2 = 0.5*m2*(v_G2.'*v_G2) + 0.5*(m2*l2*l2/12)*(w2.'*w2);
T3 = 0.5*m3*(v_G3.'*v_G3) + 0.5*(m3*l3*l3/12)*(w3.'*w3);
T4 = 0.5*m4*(v_G4.'*v_G4) + 0.5*(m4*l4*l4/12)*(w4.'*w4);
T = simplify(T1 + T2 + T3 + T4);

% Potential Energy
V = simplify(0.5*k1*(q4 - q3)^2 + 0.5*k2*(q5 - q4)^2 + 0.5*k3*(q6 - q5)^2);

% Lagrangian Function
L = simplify(T - V);

% Total Mechanical Energy
E = simplify(T + V);
% matlabFunction(E, 'File', 'Energy');

% Matrices Calculation
pl_pdq = jacobian(L,dq).';
M = jacobian(pl_pdq,dq);
M = simplify(M);

pl_pq = jacobian(L,q).';
B = diff(pl_pdq,t) + jacobian(pl_pdq,q)*dq - pl_pq;
B = simplify(B);

Q = [0;0;0;0;0;0];
F = Q - B;
F = simplify(F);

% Integrated Multipliers Method
% A_IM = [eye(6,6) zeros(6,6) zeros(6,4); zeros(6,6) M -a.'; zeros(4,6) -a zeros(4,4)];
% B_IM = [dq; F; bdot + adot*dq];
% A_IM = simplify(A_IM);
% B_IM = simplify(B_IM);
% matlabFunction(A_IM, 'File', 'A_IntMul');
% matlabFunction(B_IM, 'File', 'B_IntMul');

% Augmented Method
% A_Au = [M -a.'; -a zeros(4,4)];
% B_Au = [F; bdot + adot*dq];
% A_Au = simplify(A_Au);
% B_Au = simplify(B_Au);
% matlabFunction(A_Au, 'File', 'A_Aug');
% matlabFunction(B_Au, 'File', 'B_Aug');

% Elimination Method
% M1 = M(1:4,:);
% M2 = M(5:6,:);
% F1 = F(1:4);
% F2 = F(5:6);
% a1 = a(:,1:4);
% a2 = a(:,5:6);
% A_El = [M2 - a2.'*(inv(a1))'*M1; -a];
% B_El = [F2 - a2.'*(inv(a1))'*F1; bdot + adot*dq];
% A_El = simplify(A_El);
% B_El = simplify(B_El);
% matlabFunction(A_El, 'File', 'A_Elim');
% matlabFunction(B_El, 'File', 'B_Elim');

% Embedding Method
% dqD = [dq2; dq3; dq4; dq5];
% dqI = [dq1; dq6];     %xdot and theta4dot are chosen to be independent
% a11 = a(:,2:5);
% a22(:,1) = a(:,1);
% a22(:,2) = a(:,6);
% Wbar = [1 0; -a11\a22; 0 1];
% Xbar = [0; -a11\b; 0];
% dWbar = zeros(size(Wbar));
% dXbar = zeros(size(Xbar));
% for i = 1:length(q)
%     dWbar = dWbar + diff(Wbar, q(i))*dq(i);
%     dXbar = dXbar + diff(Xbar, q(i))*dq(i);
% end
% matlabFunction(M, 'File', 'MassMat');
% matlabFunction(F, 'File', 'QminusB');
% matlabFunction(Wbar, 'File', 'Wbar_Embd');
% matlabFunction(dWbar, 'File', 'Wbardot_Embd');
% matlabFunction(Xbar, 'File', 'Xbar_Embd');
% matlabFunction(dXbar, 'File', 'Xbardot_Embd');

%Greenwood Method
% matlabFunction(bdot, 'File', 'ConstrBiasdot');
% matlabFunction(adot, 'File', 'ConstrJacobdot');
% matlabFunction(b, 'File', 'ConstrBias');
% matlabFunction(a, 'File', 'ConstrJacob');

%Modified Lagrange Method
% e11 = [cos(q3);sin(q3);0];
% u1 = dot(v_G1,e11);
% u2 = dq3;
% u =[u1;u2];
% Y = jacobian(u,dq);
% Z = u - Y*dq;
% W = ([Y;a])\([eye(2,2);zeros(4,2)]);
% X = -([Y;a])\[Z;b];
% W = simplify(W);
% X = simplify(X);
% Wdot = zeros(size(W));
% Xdot = zeros(size(X));
% for i = 1:length(q)
%     Wdot = Wdot + diff(W, q(i))*dq(i);
%     Xdot = Xdot + diff(X, q(i))*dq(i);
% end
% Wdot = simplify(Wdot);
% Xdot = simplify(Xdot);
% Mprime = W.'*M*W;
% Bprime = W.'*M*(Wdot*u + Xdot) + W.'*B;
% Qprime = W.'*Q;
% 
% Mprime = simplify(Mprime);
% Bprime = simplify(Bprime);
% Qprime = simplify(Qprime);

% matlabFunction(Mprime, 'File', 'Mpr_Ml');
% matlabFunction(Bprime, 'File', 'Bpr_Ml');
% matlabFunction(Qprime, 'File', 'Qpr_Ml');
% matlabFunction(W, 'File', 'W_Ml');
% matlabFunction(X, 'File', 'X_Ml');
% matlabFunction(u, 'File', 'u_Ml');

% Gibss-Appel Method
e11 = [cos(q3);sin(q3);0];
u1 = dot(v_G1,e11);
u2 = dq3;
u =[u1;u2];
Y = jacobian(u,dq);
Z = u - Y*dq;
W = ([Y;a])\([eye(2,2);zeros(4,2)]);
X = -([Y;a])\[Z;b];
W = simplify(W);
X = simplify(X);

% Stating Velocites in terms of u
U=[U1;U2];     %Velocities presented here are funcitons of q and U
dq_itu = W*U+X;
v_G4_itu = subs(v_G4,{dq(1), dq(2), dq(3), dq(4), dq(5), dq(6)}, {dq_itu(1), dq_itu(2), dq_itu(3), dq_itu(4), dq_itu(5), dq_itu(6)});
v_G3_itu = subs(v_G3,{dq(1), dq(2), dq(3), dq(4), dq(5), dq(6)}, {dq_itu(1), dq_itu(2), dq_itu(3), dq_itu(4), dq_itu(5), dq_itu(6)});
v_G2_itu = subs(v_G2,{dq(1), dq(2), dq(3), dq(4), dq(5), dq(6)}, {dq_itu(1), dq_itu(2), dq_itu(3), dq_itu(4), dq_itu(5), dq_itu(6)});
v_G1_itu = subs(v_G1,{dq(1), dq(2), dq(3), dq(4), dq(5), dq(6)}, {dq_itu(1), dq_itu(2), dq_itu(3), dq_itu(4), dq_itu(5), dq_itu(6)});
w1_itu = subs(w1,{dq(1), dq(2), dq(3), dq(4), dq(5), dq(6)}, {dq_itu(1), dq_itu(2), dq_itu(3), dq_itu(4), dq_itu(5), dq_itu(6)});
w2_itu = subs(w2,{dq(1), dq(2), dq(3), dq(4), dq(5), dq(6)}, {dq_itu(1), dq_itu(2), dq_itu(3), dq_itu(4), dq_itu(5), dq_itu(6)});
w3_itu = subs(w3,{dq(1), dq(2), dq(3), dq(4), dq(5), dq(6)}, {dq_itu(1), dq_itu(2), dq_itu(3), dq_itu(4), dq_itu(5), dq_itu(6)});
w4_itu = subs(w4,{dq(1), dq(2), dq(3), dq(4), dq(5), dq(6)}, {dq_itu(1), dq_itu(2), dq_itu(3), dq_itu(4), dq_itu(5), dq_itu(6)});

v_G4_itu = simplify(v_G4_itu);     
v_G3_itu = simplify(v_G3_itu);
v_G2_itu = simplify(v_G2_itu);
v_G1_itu = simplify(v_G1_itu);
w1_itu = simplify(w1_itu);
w2_itu = simplify(w2_itu);
w3_itu = simplify(w3_itu);
w4_itu = simplify(w4_itu);

%Calculating Angular Accelerations 
dU = [dU1;dU2];

alpha1 = zeros(3,1);
for i=1:length(q)
alpha1 = alpha1 + diff(w1_itu,q(i))*dq_itu(i);
end
alpha1 = alpha1 + diff (w1_itu, U1)*dU1 + diff (w1_itu, U2)*dU2;
alpha1 = simplify(alpha1);

alpha2 = zeros(3,1);
for i=1:length(q)
alpha2 = alpha2 + diff(w2_itu,q(i))*dq_itu(i);
end
alpha2 = alpha2 + diff (w2_itu, U1)*dU1 + diff (w2_itu, U2)*dU2;
alpha2 = simplify(alpha2);

alpha3 = zeros(3,1);
for i=1:length(q)
alpha3 = alpha3 + diff(w3_itu,q(i))*dq_itu(i);
end
alpha3 = alpha3 + diff (w3_itu, U1)*dU1 + diff (w3_itu, U2)*dU2;
alpha3 = simplify(alpha3);

alpha4 = zeros(3,1);
for i=1:length(q)
alpha4 = alpha4 + diff(w4_itu,q(i))*dq_itu(i);
end
alpha4 = alpha4 + diff (w4_itu, U1)*dU1 + diff (w4_itu, U2)*dU2;
alpha4 = simplify(alpha4);

%Calculating Linear Accelerations
accel1 = zeros(3,1);
for i=1:length(q)
accel1 = accel1+ diff(v_G1_itu,q(i))*dq_itu(i);
end
accel1 = accel1+ diff(v_G1_itu,U1)*dU1 + diff(v_G1_itu,U2)*dU2;
accel1 = simplify(accel1);

accel2 = zeros(3,1);
for i=1:length(q)
accel2 = accel2+ diff(v_G2_itu,q(i))*dq_itu(i);
end
accel2 = accel2+ diff(v_G2_itu,U1)*dU1 + diff(v_G2_itu,U2)*dU2;
accel2 = simplify(accel2);

accel3 = zeros(3,1);
for i=1:length(q)
accel3 = accel3+ diff(v_G3_itu,q(i))*dq_itu(i);
end
accel3 = accel3+ diff(v_G3_itu,U1)*dU1 + diff(v_G3_itu,U2)*dU2;
accel3 = simplify(accel3);

accel4 = zeros(3,1);
for i=1:length(q)
accel4 = accel4+ diff(v_G4_itu,q(i))*dq_itu(i);
end
accel4 = accel4+ diff(v_G4_itu,U1)*dU1 + diff(v_G4_itu,U2)*dU2;
accel4 = simplify(accel4);

% Partial Velocites
v_G1_1 = diff(v_G1_itu, U1);
v_G1_2 = diff(v_G1_itu, U2);
w1_1 = diff(w1_itu, U1);
w1_2 = diff(w1_itu, U2);

v_G2_1 = diff(v_G2_itu, U1);
v_G2_2 = diff(v_G2_itu, U2);
w2_1 = diff(w2_itu, U1);
w2_2 = diff(w2_itu, U2);

v_G3_1 = diff(v_G3_itu, U1);
v_G3_2 = diff(v_G3_itu, U2);
w3_1 = diff(w3_itu, U1);
w3_2 = diff(w3_itu, U2);

v_G4_1 = diff(v_G4_itu, U1);
v_G4_2 = diff(v_G4_itu, U2);
w4_1 = diff(w4_itu, U1);
w4_2 = diff(w4_itu, U2);

% Generalized Forces in terms of u
M1 = -k1*(q4 - q3)*[0;0;1];
M2 = -k2*(q5 - q4)*[0;0;1];
M3 = -k3*(q6 - q5)*[0;0;1];

Qpr1 = dot(M1,w2_1) + dot(-M1,w1_1) + dot(M2,w3_1) + dot(-M2,w2_1) + dot(M3,w4_1) + dot(-M3,w3_1);
Qpr2 = dot(M1,w2_2) + dot(-M1,w1_2) + dot(M2,w3_2) + dot(-M2,w2_2) + dot(M3,w4_2) + dot(-M3,w3_2);
Qpr = [Qpr1;Qpr2];

% Building Gibss Function
I1 = 1/12*m1*l1*l1; I2 = 1/12*m2*l2*l2; I3 = 1/12*m3*l3*l3; I4 = 1/12*m4*l4*l4;
S = 0.5*(m1*(dot(accel1,accel1)) + m2*(dot(accel2,accel2)) + m3*(dot(accel3,accel3)) + m4*(dot(accel4,accel4))) ...
    + 0.5*(dot(alpha1,I1*alpha1) + (dot(alpha2,I2*alpha2)) + (dot(alpha3,I3*alpha3)) + (dot(alpha4,I4*alpha4))) ...
    + dot(alpha1,(cross(w1_itu, (I1*w1_itu)))) + dot(alpha2,(cross(w2_itu, (I2*w2_itu)))) + dot(alpha3,(cross(w3_itu, (I3*w3_itu)))) + dot(alpha4,(cross(w4_itu, (I4*w4_itu))));

% Equations of Motion
Sdot1 = diff(S,dU1);
Sdot2 = diff(S,dU2);
eq = [Sdot1;Sdot2] - Qpr;
M_g = jacobian(eq,dU);
B_g = eq - M_g*dU;
% matlabFunction(M_g, 'File', 'M_gibss');
% matlabFunction(B_g, 'File', 'B_gibss');

% Kane Method
dH1 = I1*alpha1 + cross(w1_itu, I1*w1_itu);
dH2 = I2*alpha2 + cross(w2_itu, I2*w2_itu);
dH3 = I3*alpha3 + cross(w3_itu, I3*w3_itu);
dH4 = I4*alpha4 + cross(w4_itu, I4*w4_itu);

% Building Fstars
Fstar1 = -(m1*dot(accel1,v_G1_1)+m2*dot(accel2,v_G2_1)+m3*dot(accel3,v_G3_1)+m4*dot(accel4,v_G4_1))...
    - (dot(dH1,w1_1)+dot(dH2,w2_1)+dot(dH3,w3_1)+dot(dH4,w4_1));

Fstar2 = -(m1*dot(accel1,v_G1_2)+m2*dot(accel2,v_G2_2)+m3*dot(accel3,v_G3_2)+m4*dot(accel4,v_G4_2))...
    - (dot(dH1,w1_2)+dot(dH2,w2_2)+dot(dH3,w3_2)+dot(dH4,w4_2));

% Equations of Motion
eq1_kane = Fstar1+Qpr1;
eq2_kane = Fstar2+Qpr2;

eq_kane = [eq1_kane;eq2_kane];

M_k = jacobian(eq_kane,dU);
B_k = eq_kane - M_k*dU;
% matlabFunction(M_k, 'File', 'M_kane');
% matlabFunction(B_k, 'File', 'B_kane');

