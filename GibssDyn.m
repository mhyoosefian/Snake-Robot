function dz = GibssDyn(t,z)
dz = z;
q1 = z(1) ; q2 = z(2); q3 = z(3) ; q4 = z(4); q5 = z(5) ; q6 = z(6); 
U1 = z(7); U2 = z(8);
U = [U1;U2];
dU1 = 5; dU2 = 0; 
global m1 l1 m2 l2 m3 l3 m4 l4 k1 k2 k3
M_g = M_gibss(l1,m1,m2,m3,m4,q3,q4,q5,q6);
B_g = B_gibss(U1,U2,dU1,dU2,k1,k2,k3,l1,l2,l3,l4,m1,m2,m3,m4,q3,q4,q5,q6);

du = -inv(M_g)*B_g;
W = W_Ml(l1,l2,l3,l4,q3,q4,q5,q6);
X = X_Ml;
dq = W*U+X;

dz=[dq;du];