function D=HJB_Function(~,Y)
%% Variable definition

global J V X0 D1 D2 dLq A JacobA
ST=1e-3;  Par_D1=D1;Par_D2=D2;
Par_dLq=dLq;Par_A=A; Par_JacobA=JacobA;

%%  Variable Calculations

s1=Y(1);s2=Y(2);s3=Y(3);ds1=Y(4);ds2=Y(5);ds3=Y(6);
L1=Y(7);L2=Y(8);L3=Y(9);L4=Y(10);L5=Y(11);L6=Y(12);
S=[s1;s2;s3];
options=optimset('Display','off','TolFun',1e-10);
X0=fsolve(@kinematic_Function,X0,options,S);
x=X0(1);y=X0(2);z=X0(3);say=X0(4);tet=X0(5);
fii=X0(6);alfa1=X0(7);alfa2=X0(8);alfa3=X0(9);
A0_OR=eval(A);S_omit_OR=null(A0_OR,'r');
V=[ds1;ds2;ds3];  dq_OR=S_omit_OR*V;
dx=dq_OR(1);dy=dq_OR(2);dz=dq_OR(3);
dsay=dq_OR(4);dtet=dq_OR(5);dfii=dq_OR(6);
dalfa1=dq_OR(7);dalfa2=dq_OR(8);dalfa3=dq_OR(9);
Variables=[s1 s2 s3 ds1 ds2 ds3]';   Par_X0=X0;

%%  Original Dinamics Equations

JacobA0_OR=eval(JacobA); DA=reshape(JacobA0_OR*dq_OR,9,12);
DS_omit_OR=-pinv(A0_OR)*DA*S_omit_OR;
M_OR=eval(D2);C=eval(D1);G_OR=eval(dLq);
L=[L1;L2;L3;L4;L5;L6];X=[s1;s2;s3;ds1;ds2;ds3];
InvM0_OR=inv(S_omit_OR'*M_OR*S_omit_OR);
F_OR=-InvM0_OR'*L(4:6);
h=InvM0_OR*S_omit_OR'*(G_OR-C*dq_OR-M_OR*DS_omit_OR*V);
Deq=[V;h]+[zeros(3);InvM0_OR]*F_OR;

%%  Equations

Step=[-2*ST; -ST;  ST;  2*ST];

for i=1:6
Dif=Variables(i)*ones(4,1)+Step;   

if 3<i
%------------------------------------------------------------
Y0=[Variables(1:i-1);Dif(1);Variables(i+1:6)];
options=optimset('Display','off','TolFun',1e-10);
X1=fsolve(@kinematic_Function,Par_X0,options,Y0(1:3));
H1=Parallel_Function(Y0,X1,F_OR,L,Par_D1,Par_D2,Par_dLq,Par_A,Par_JacobA);

Y0=[Variables(1:i-1);Dif(2);Variables(i+1:6)];
H2=Parallel_Function(Y0,X1,F_OR,L,Par_D1,Par_D2,Par_dLq,Par_A,Par_JacobA);

Y0=[Variables(1:i-1);Dif(3);Variables(i+1:6)];
H3=Parallel_Function(Y0,X1,F_OR,L,Par_D1,Par_D2,Par_dLq,Par_A,Par_JacobA);

Y0=[Variables(1:i-1);Dif(4);Variables(i+1:6)];
H4=Parallel_Function(Y0,X1,F_OR,L,Par_D1,Par_D2,Par_dLq,Par_A,Par_JacobA);
%------------------------------------------------------------

else
    
%------------------------------------------------------------
Y0=[Variables(1:i-1);Dif(1);Variables(i+1:6)];
options=optimset('Display','off','TolFun',1e-10);
X1=fsolve(@kinematic_Function,Par_X0,options,Y0(1:3));
H1=Parallel_Function(Y0,X1,F_OR,L,Par_D1,Par_D2,Par_dLq,Par_A,Par_JacobA);

Y0=[Variables(1:i-1);Dif(2);Variables(i+1:6)];
options=optimset('Display','off','TolFun',1e-10);
X1=fsolve(@kinematic_Function,Par_X0,options,Y0(1:3));
H2=Parallel_Function(Y0,X1,F_OR,L,Par_D1,Par_D2,Par_dLq,Par_A,Par_JacobA);

Y0=[Variables(1:i-1);Dif(3);Variables(i+1:6)];
options=optimset('Display','off','TolFun',1e-10);
X1=fsolve(@kinematic_Function,Par_X0,options,Y0(1:3));
H3=Parallel_Function(Y0,X1,F_OR,L,Par_D1,Par_D2,Par_dLq,Par_A,Par_JacobA);

Y0=[Variables(1:i-1);Dif(4);Variables(i+1:6)];
options=optimset('Display','off','TolFun',1e-10);
X1=fsolve(@kinematic_Function,Par_X0,options,Y0(1:3));
H4=Parallel_Function(Y0,X1,F_OR,L,Par_D1,Par_D2,Par_dLq,Par_A,Par_JacobA);
%------------------------------------------------------------
end

DiFF_X(i,1)=(-H4+8*H3-8*H2+H1)/(12*ST);
end

D=[Deq;-DiFF_X];    J=J+1; disp(J);
end
