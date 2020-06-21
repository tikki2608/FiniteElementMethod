syms r0a r1a ka Ia La rhoa xa ra c1 c2
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
% solving for integration constant
c1=solve(dTdxa(r0,r1,k,I,rho,L,x_dTdx,c1));
c2=solve(Ta(r0,r1,k,I,rho,L,x_T,c1,c2)==20,c2);
% Temperatre and Heat Flux calcuation
x=linspace(0,0.01,100);
Temp=Ta(r0,r1,k,I,rho,L,x,c1,c2);
Heat_flux=-HF(r0,r1,k,I,rho,L,x,c1);
% Plot of Temperature and Heat Flux 
figure('name',"Temprature Profile");
plot(x,Temp)
xlim([0 L])
title('Temperature Profile')
xlabel('x distance (m)')
ylabel('Temeperature (^0C)')
figure('name',"Heat Flux Profile");
plot(x,Heat_flux)
xlim([0 L])
title('Heat Flux Profile')
xlabel('x distance (m)')
ylabel('Heat Flux (W/m^2)')
