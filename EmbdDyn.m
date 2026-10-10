function dz = EmbdDyn(t,z)
dz = z;
global m1 l1 m2 l2 m3 l3 m4 l4 k1 k2 k3

q1 = z(1); q2 = z(2); q3 = z(3); q4 = z(4); q5 = z(5); q6 = z(6);
dq1 = z(7); dq6 = z(8);

dqI = [dq1;dq6];

Wbar = Wbar_Embd(l1,l2,l3,l4,q3,q4,q5,q6);
Xbar = Xbar_Embd;

qdot = Wbar*dqI + Xbar;
dq2 = qdot(2); dq3 = qdot(3); dq4 = qdot(4); dq5 = qdot(5);

dWbar = Wbardot_Embd(dq3,dq4,dq5,dq6,l1,l2,l3,l4,q3,q4,q5,q6);
dXbar = Xbardot_Embd;

M = MassMat(l1,l2,l3,l4,m1,m2,m3,m4,q3,q4,q5,q6);
F = QminusB(dq3,dq4,dq5,dq6,k1,k2,k3,l1,l2,l3,l4,m2,m3,m4,q3,q4,q5,q6);


A_Em = Wbar.'*M*Wbar;
B_Em = Wbar.'*(F - M*(dWbar*dqI+dXbar));

dz = [qdot; A_Em\B_Em];