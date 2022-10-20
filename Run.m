clear all;clc

%% Definitions

global A D2 X0 J
T0=0; Tf=1; Mesh=5;
J=0;   S=[0.4;0.4;0.4];   V=[0;0;0];
X0=[0;0.5;0;0;0;0;1.2;1.2;1.2];

%% Calling functions

Parametric_Calculations;

%% Solve

options=optimset('Display','off','TolFun',1e-5);
X0=fsolve(@kinematic_Function,X0,options,S);
X0_OR=X0;
options = bvpset('RelTol',1e-3,'AbsTol',1e-4);
solinit = bvpinit(linspace(T0,Tf,Mesh),@Initial_Conditions);
sol = bvp5c(@HJB_Function,@Boundary_Conditions,solinit,options);
xint = linspace(T0,Tf,Mesh);
Sxint = deval(sol,xint);

%%  Calculating Forces

for i=1:Mesh
S=Sxint(1:3,i); L=Sxint(10:12,i);
options=optimset('Display','off','TolFun',1e-10);
X0_OR=fsolve(@kinematic_Function,X0_OR,options,S);
x=X0_OR(1);y=X0_OR(2);z=X0_OR(3);say=X0_OR(4);tet=X0_OR(5);fii=X0_OR(6);
alfa1=X0_OR(7);alfa2=X0_OR(8);alfa3=X0_OR(9);
A0=eval(A);S_omit=null(A0,'r');
M=eval(D2);  InvM0=inv(S_omit'*M*S_omit);
Coordinates(:,i)=X0_OR;
F(:,i)=-InvM0'*L;
end

%% Plots

plot(xint,Sxint(1:3,:),'linewidth',2);legend('S^*_1','S^*_2','S^*_3');xlabel('Time(s)','fontsize',14);ylabel('Position(m)','fontsize',14);grid;
figure
plot(xint,Sxint(4:6,:),'linewidth',2);legend('dS^*_1','dS^*_2','dS^*_3');xlabel('Time(s)','fontsize',14);ylabel('Spees(m/s)','fontsize',14);grid;
figure
plot(xint,F,'linewidth',2);xlabel('Time(s)','fontsize',14);ylabel('Forces(N)','fontsize',14);legend('F^*_1','F^*_2','F^*_3','location','best');grid on;

