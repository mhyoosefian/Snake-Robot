function dz = ModlagDyn(t,z)
dz = z;
global m1 l1 m2 l2 m3 l3 m4 l4 k1 k2 k3

q1 = z(1); q2 = z(2); q3 = z(3); q4 = z(4); q5 = z(5); q6 = z(6);
u1 = z(7); u2 = z(8);
u = [u1;u2];
W = W_Ml(l1,l2,l3,l4,q3,q4,q5,q6);
X = X_Ml;
qdot = W*u + X;
dq1 = qdot(1); dq2 = qdot(2); dq3 = qdot(3); dq4 = qdot(4); dq5 = qdot(5); dq6 = qdot(6);
Qprime = Qpr_Ml;
Mprime = Mpr_Ml(l1,m1,m2,m3,m4,q3,q4,q5,q6); 
Bprime = Bpr_Ml(dq1,dq2,dq3,dq4,dq5,dq6,k1,k2,k3,l1,l2,l3,l4,m1,m2,m3,m4,q3,q4,q5,q6);

A = W*u + X;
B = Mprime\(Qprime - Bprime);


dz = [A;B];