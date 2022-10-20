%%  Variable Definition

a=0.8;b=0.2;l=0.5;g=9.8;M=1;m=0.1;
syms x y z say tet fii alfa1 alfa2 alfa3 s1 s2 s3
syms dx dy dz dsay dtet dfii dalfa1 dalfa2 dalfa3 ds1 ds2 ds3
syms d2x d2y d2z d2say d2tet d2fii d2alfa1 d2alfa2 d2alfa3 d2s1 d2s2 d2s3
syms T1 T2 T3
global Constraint A JacobA D1 D2  dLq

%%  Kinematic

o=[0 0 0]';
a1=[0 0 a]';
a2=[a*sqrt(3)/2 0 -a/2]';
a3=[-a*sqrt(3)/2 0 -a/2]';
Hc1=[0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 1 0 0;0 0 1 a-s1;0 0 0 1]*...
    [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1];
Hc2=[cosd(120) 0 sind(120) 0;0 1 0 0;-sind(120) 0 cosd(120) 0;0 0 0 1]*...
    [0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 1 0 0;0 0 1 a-s2;0 0 0 1]*...
    [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1];
Hc3=[cosd(-120) 0 sind(-120) 0;0 1 0 0;-sind(-120) 0 cosd(-120) 0;0 0 0 1]*...
    [0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 1 0 0;0 0 1 a-s3;0 0 0 1]*...
    [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1];

c1=Hc1*[0;0;0;1];
c2=Hc2*[0;0;0;1];
c3=Hc3*[0;0;0;1];
b1=Hc1*[l*sin(alfa1);l*cos(alfa1);0;1];
b2=Hc2*[l*sin(alfa2);l*cos(alfa2);0;1];
b3=Hc3*[l*sin(alfa3);l*cos(alfa3);0;1];
b1=b1(1:3);b2=b2(1:3);b3=b3(1:3);
Rsay=[1 0 0;0 cos(say) -sin(say);0 sin(say) cos(say)];
Rfi=[cos(fii) -sin(fii) 0;sin(fii) cos(fii) 0;0 0 1];
Rtet=[cos(tet) 0 sin(tet);0 1 0;-sin(tet) 0 cos(tet)];
R=Rfi*Rtet*Rsay;
B1=[x;y;z]+R*[0;0;b];
B2=[x;y;z]+R*[b*sqrt(3)/2;0;-b/2];
B3=[x;y;z]+R*[-b*sqrt(3)/2;0;-b/2];

%%  Kinematic Constraint Equations

f1=b1(1)-B1(1);
f2=b1(2)-B1(2);
f3=b1(3)-B1(3);
f4=b2(1)-B2(1);
f5=b2(2)-B2(2);
f6=b2(3)-B2(3);
f7=b3(1)-B3(1);
f8=b3(2)-B3(2);
f9=b3(3)-B3(3);
Constraint=[f1;f2;f3;f4;f5;f6;f7;f8;f9];
A=jacobian(Constraint,[x y z say tet fii alfa1 alfa2 alfa3 s1 s2 s3]);
JacobA=jacobian(reshape(A,9*12,1),[x y z say tet fii alfa1 alfa2 alfa3 s1 s2 s3]);

%%  Dynamics Equations

Ip=R*[1/4*M*b^2,0,0;0,1/2*M*b^2,0;0,0,1/4*M*b^2]*R.';
IL1=Hc1(1:3,1:3)*[cos(-alfa1) -sin(-alfa1) 0;sin(-alfa1) cos(-alfa1) 0;0 0 1]*...
    [1/3*m*l^2,0,0;0,0,0;0,0,1/3*m*l^2]*[cos(-alfa1) -sin(-alfa1) 0;sin(-alfa1) cos(-alfa1) 0;0 0 1].'*Hc1(1:3,1:3).';
IL2=Hc2(1:3,1:3)*[cos(-alfa2) -sin(-alfa2) 0;sin(-alfa2) cos(-alfa2) 0;0 0 1]*...
    [1/3*m*l^2,0,0;0,0,0;0,0,1/3*m*l^2]*[cos(-alfa2) -sin(-alfa2) 0;sin(-alfa2) cos(-alfa2) 0;0 0 1].'*Hc2(1:3,1:3).';
IL3=Hc3(1:3,1:3)*[cos(-alfa3) -sin(-alfa3) 0;sin(-alfa3) cos(-alfa3) 0;0 0 1]*[1/3*m*l^2,0,0;0,0,0;0,0,1/3*m*l^2]*...
    [cos(-alfa3) -sin(-alfa3) 0;sin(-alfa3) cos(-alfa3) 0;0 0 1].'*Hc3(1:3,1:3).';
omegaL1=Hc1(1:3,1:3)*[0;0;-dalfa1];
omegaL2=Hc2(1:3,1:3)*[0;0;-dalfa2];
omegaL3=Hc3(1:3,1:3)*[0;0;-dalfa3];
omegap=[dsay;dtet;dfii];
Vp=[dx;dy;dz];
G1=Hc1*[0.5*l*sin(alfa1);0.5*l*cos(alfa1);0;1];
G2=Hc2*[0.5*l*sin(alfa2);0.5*l*cos(alfa2);0;1];
G3=Hc3*[0.5*l*sin(alfa3);0.5*l*cos(alfa3);0;1];
VL1=jacobian(G1(1:3),[s1 alfa1])*[ds1;dalfa1];
VL2=jacobian(G2(1:3),[s2 alfa2])*[ds2;dalfa2];
VL3=jacobian(G3(1:3),[s3 alfa3])*[ds3;dalfa3];
Tp=1/2*omegap.'*Ip*omegap+1/2*M*(Vp.'*Vp);
TL1=1/2*omegaL1.'*IL1*omegaL1+1/2*m*(VL1.'*VL1);
TL2=1/2*omegaL2.'*IL2*omegaL2+1/2*m*(VL2.'*VL2);
TL3=1/2*omegaL3.'*IL3*omegaL3+1/2*m*(VL3.'*VL3);
TL=TL1+TL2+TL3;
U=M*g*y+m*g*l/2*(sin(alfa1)+sin(alfa2)+sin(alfa3));
L=Tp+TL-U;
dLq=jacobian(L,[x y z say tet fii alfa1 alfa2 alfa3 s1 s2 s3]).';
dLdq=jacobian(L,[dx dy dz dsay dtet dfii dalfa1 dalfa2 dalfa3 ds1 ds2 ds3]).';
D2=jacobian(dLdq,[dx dy dz dsay dtet dfii dalfa1 dalfa2 dalfa3 ds1 ds2 ds3]);
D1=jacobian(dLdq,[x y z say tet fii alfa1 alfa2 alfa3 s1 s2 s3]);


A=simplify(A);
JacobA=simplify(JacobA);
dLq=simplify(dLq);
D2=simplify(D2);
D1=simplify(D1);
