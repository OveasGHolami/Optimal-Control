function [out]=Parallel_Function(Y0,X1,F_OR,L,D1,D2,dLq,A,JacobA)
%% Variable definition

R=eye(6);P=eye(3);
s1=Y0(1);s2=Y0(2);s3=Y0(3); ds1=Y0(4);ds2=Y0(5);ds3=Y0(6);
x=X1(1);y=X1(2);z=X1(3);say=X1(4);tet=X1(5);
fii=X1(6);alfa1=X1(7);alfa2=X1(8);alfa3=X1(9);
A0=eval(A);  S_omit=null(A0,'r');
V=[ds1;ds2;ds3]; dq=S_omit*V;
dx=dq(1);dy=dq(2);dz=dq(3);dsay=dq(4);dtet=dq(5);
dfii=dq(6);dalfa1=dq(7);dalfa2=dq(8);dalfa3=dq(9);

%%  Equations

JacobA0=eval(JacobA);  DA=reshape(JacobA0*dq,9,12);
DS_omit=-pinv(A0)*DA*S_omit;
M=eval(D2);G=eval(dLq);C=eval(D1);
X=[s1;s2;s3;ds1;ds2;ds3];
InvM0=inv(S_omit'*M*S_omit);
h=InvM0*S_omit'*(G-C*dq-M*DS_omit*V);
out=0.5*(X'*R*X+F_OR'*P*F_OR)+...
    L'*([V;h]+[zeros(3);InvM0]*F_OR);
end
