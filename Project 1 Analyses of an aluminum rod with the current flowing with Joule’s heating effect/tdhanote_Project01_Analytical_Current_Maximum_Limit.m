syms r0a r1a ka Ia La rhoa xa ra c1 c2 C1 C2 
% parameters equations
ra=r1a +(r0a-r1a)*(xa/La);
Aa=pi()*ra^2;
%Derivation
pa=int(1/Aa,xa);
HF(r0a,r1a,ka,Ia,rhoa,La,xa,c1)=-(pa*Ia^2*rhoa)/Aa-c1/Aa;
dTdxa(r0a,r1a,ka,Ia,rhoa,La,xa,c1)=-((pa*Ia^2*rhoa)/(ka*Aa))-c1/(ka*Aa);
Ta(r0a,r1a,ka,Ia,rhoa,La,xa,c1,c2)=int(dTdxa,xa)+c2;
%Inputs
k=205;
r0=0.002;
r1=0.001;
L=0.01;
I=1000;
rho=2.82e-8;
x_T=L;
x_dTdx=0;
c1_I_max=solve(dTdxa(r0,r1,k,Ia,rho,L,x_dTdx,C1),C1);
c2_I_max=solve(Ta(r0,r1,k,Ia,rho,L,x_T,c1_I_max,C2)==20,C2);
I_max=solve(Ta(r0, r1, k, Ia, rho, L, xp, c1_I_max, c2_I_max)==660,Ia);
I_max=max(double(I_max));
disp(['Maximum Current,Imax = ',num2str(I_max),' A'])