function dz = GreenDyn(t,z)

dz = z;
global m1 l1 m2 l2 m3 l3 m4 l4 k1 k2 k3

q1 = z(1); q2 = z(2); q3 = z(3); q4 = z(4); q5 = z(5); q6 = z(6);
dq1 = z(7); dq2 = z(8); dq3 = z(9); dq4 = z(10); dq5 = z(11); dq6 = z(12);

M = MassMat(l1,l2,l3,l4,m1,m2,m3,m4,q3,q4,q5,q6);
F = QminusB(dq3,dq4,dq5,dq6,k1,k2,k3,l1,l2,l3,l4,m2,m3,m4,q3,q4,q5,q6);
a = ConstrJacob(l1,l2,l3,l4,q3,q4,q5,q6);
b = ConstrBias;
adot = ConstrJacobdot(dq3,dq4,dq5,dq6,l1,l2,l3,q3,q4,q5,q6);
bdot = ConstrBiasdot;

A = (a.'/(a*(M\(a.'))));
qdot = [dq1;dq2;dq3;dq4;dq5;dq6];
B = -bdot - adot*qdot - a*(M\F);
C = M\(F + A*B);
dz = [qdot;C];