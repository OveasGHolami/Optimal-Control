function [C]=Constraint(Coordinats)
%%  Variable Definition

a=0.8;b=0.2;l=0.5;
x=Coordinats(1); y=Coordinats(2); z=Coordinats(3); say=Coordinats(4);
tet=Coordinats(5); fii=Coordinats(6); alfa1=Coordinats(7); alfa2=Coordinats(8);
alfa3=Coordinats(9); s1=Coordinats(10); s2=Coordinats(11); s3=Coordinats(12);

%%  Kinematic

Hc1=[0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 1 0 0;0 0 1 a-s1;0 0 0 1]*...
    [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1];
Hc2=[cosd(120) 0 sind(120) 0;0 1 0 0;-sind(120) 0 cosd(120) 0;0 0 0 1]*...
    [0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 1 0 0;0 0 1 a-s2;0 0 0 1]*...
    [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1];
Hc3=[cosd(-120) 0 sind(-120) 0;0 1 0 0;-sind(-120) 0 cosd(-120) 0;0 0 0 1]*...
    [0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 1 0 0;0 0 1 a-s3;0 0 0 1]*...
    [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1];

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
C=[f1;f2;f3;f4;f5;f6;f7;f8;f9];
end
