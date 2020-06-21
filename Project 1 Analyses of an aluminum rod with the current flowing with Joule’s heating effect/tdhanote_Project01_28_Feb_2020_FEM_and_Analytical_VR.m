function [fe,K,T,I_max] = tdhanote_Project01_28_Feb_2020_FEM_and_Analytical_VR(nn_per_element,number_of_elements,k,r0, r1,L,T0,no_qpts)
%% input defination 

% nn_per_element = Number of nodes per element
% number_of_elements = Number of elements
% k = Thermal conductivity
% r0 = Radius at the ends of bar
% r1 = radius at the center of bar
% roh_0 = Electrical Resistivity is constant or vaariable
% no_qpts = Number of quadrature Points

%% Inpts
if nargin==0
nn_per_element=2;  
number_of_elements=30;
k=205;
r0=0.002;
r1=0.001;
L=0.01;
I=1000;
no_qpts=10;
T0=20;
end

%% Heat Generation as a function of x

s=@(x,rho)(((I^2)*rho)/(pi()*(r1 + (r0-r1)*(x/L))^2));

%% Arear as a function of x

A=@(x)(pi()*(r1 + (r0-r1)*(x/L))^2);

%% Quadture points selection

qpts=quadrature(no_qpts);

%% Total nodes and connectivity matrix determination Corresponding shape function

if nn_per_element==2
    nodes_total=number_of_elements+1;
    conn=[1:number_of_elements;2:number_of_elements+1];
    shape=@shape2;
elseif nn_per_element==3
     nodes_total=2*number_of_elements+1;
     conn=[1:2:2*number_of_elements-1;2:2:2*number_of_elements;3:2:2*number_of_elements+1];
     shape=@shape3;
else
 conn=[1:3:3*number_of_elements-2;2:3:3*number_of_elements-1;3:3:3*number_of_elements;4:3:3*number_of_elements+1];
 nodes_total=3*number_of_elements+1;
 shape=@shape4;
end
%% Electrical Resistivity
rho_fun=@(T_in)(2.6e-8+1.1e-10*T_in);
%% Coordiantes of each Nodes
x=linspace(0, L,nodes_total);

%% Stiffness matrix  K Calculation

K=zeros(nodes_total);
for c=conn
    xe=x(:,c);
    Ke=zeros(length(c));
    for q=qpts
        [N,dNdp]=shape(q(1));
        J=xe*dNdp;
        B=dNdp/J;
        xq=xe*N;
        Ke=Ke+B*k*A(xq)*B'*J*q(2);
    end
    sctr=c;
    K(sctr,sctr)=K(sctr,sctr)+Ke;
end

%% Heat Flux Matrix calculations

T1=ones(nodes_total,1)*T0;
T_0=ones(nodes_total,1)*T0;
Rho=rho_fun(T_0);
dT=100;
count=0;
while max(dT) >= 0.00001
    fe=zeros(nodes_total,1);
for c=conn
    xe=x(:,c);
    RHO=Rho(c);
        for i=1:no_qpts   
    [N,DNDp]=shape(qpts(1,i));
    J=xe*DNDp;
    w=qpts(2,i);
    xp=xe*N;
    rho=RHO'*N;
    F=s(xp,rho)*N*J*w;
    fe(c)=fe(c)+F;
        end
end
%% Temperature Calculation at nodes
K(nodes_total, :)= 0;
K(nodes_total,nodes_total)=1;
fe(nodes_total,1)=T0;
T=K\fe;
Rho=rho_fun(T);
dT=abs(T-T1);
T1=T;
count=count+1;
%% convergence of Temperature for Electrical resistivity change
plot(count,max(dT./T),'o')
hold on
end

%% Heat Flux do (it in terms of Paraent cooridinates)

Heat_flux=zeros(nodes_total,1);
for j=2:nodes_total
    dx=x(j-1)-x(j);
    Heat_flux(j)=k*(T(j)-T(j-1))/dx;
end
%% Analyical Soultion

syms r0a r1a ka Ia La rhoa xa ra c1 c2 C1 C2 
ra=r1a +(r0a-r1a)*(xa/La);
Aa=pi()*ra^2;
pa=int(1/Aa,xa);
HF(r0a,r1a,ka,Ia,rhoa,La,xa,c1)=-(pa*Ia^2*rhoa/Aa)-(c1/Aa);
dTdxa(r0a,r1a,ka,Ia,rhoa,La,xa,c1)=-((pa*Ia^2*rhoa)/(ka*Aa))-c1/(ka*Aa);
Ta(r0a,r1a,ka,Ia,rhoa,La,xa,c1,c2)=int(dTdxa,xa)+c2;
x_T=L;
x_dTdx=0;
c1_I_max=solve(dTdxa(r0,r1,k,Ia,rho,L,x_dTdx,C1),C1);
c2_I_max=solve(Ta(r0,r1,k,Ia,rho,L,x_T,c1_I_max,C2)==T0,C2);
I_max=solve(Ta(r0, r1, k, Ia, rho, L, x_dTdx, c1_I_max, c2_I_max)==660,Ia);
I_max=max(double(I_max));
c1=solve(dTdxa(r0,r1,k,I,rho,L,x_dTdx,c1));
c2=solve(Ta(r0,r1,k,I,rho,L,x_T,c1,c2)==T0,c2);
c1=double(c1);
c2=double(c2);
x_a=linspace(0,0.01,100);
Temp_analytical=Ta(r0,r1,k,I,rho,L,x_a,c1,c2);
Heat_flux_analytical=-HF(r0,r1,k,I,rho,L,x_a,c1);

%% Plots

figure('name',"Temprature Profile");
plot(x,T,'o--')
hold on
plot(x_a,Temp_analytical)
xlim([0 L])
title('Temperature Profile')
xlabel('X Distance')
ylabel('Temeperature')
legend('FEM','Analytical')
hold off
figure('name',"Heat Flux Profile");
plot(x,Heat_flux,'o--')
hold on
plot(x_a,Heat_flux_analytical)
xlim([0 L])
title('Heat Flux Profile')
xlabel('X distance')
ylabel('Heat Flux')
legend('FEM','Analytical')
hold off

%% Convergence

end

%% Shape functions for different elements
function [N,dNdp] = shape2(p)
N=0.5*[1-p(1);1+p(1)];
dNdp=[-0.5;0.5];
end
function [N,dNdp] = shape3(p)
N=[p(1)*(1-p(1))*-0.5;1-p(1)^2;0.5*p(1)*(p(1)+1)];
dNdp = [(p(1)-0.5);(-2*p(1));(p(1)+0.5)];
end
function [N,dNdp] = shape4(p)
N=[(p(1)^2-1/9)*(p(1)-1)*(-9/16);(p(1)^2-1)*(p(1)-1/3)*(27/16);(p(1)^2-1)*(p(1)+1/3)*(-27/16);(p(1)^2-1/9)*(p(1)+1)*(9/16)];
dNdp=[(-9/16)*(3*p(1)^2-2*p(1)-1/9);(27/16)*(3*p(1)^2-(2/3)*p(1)-1);(-27/16)*(3*p(1)^2+(2/3)*p(1)-1);(9/16)*(3*p(1)^2+2*p(1)-1/9)];
end
%% Gaussian Quadrature 
function [qpts] = quadrature(n)
%   QUADRATURE
%     quadrature(n) returns a quadrature table for a rule with n
%     integration points.  The first row of the table gives the quadrature
%     point location, and the second gives the quadrature weights.

    u = 1:n-1;
    u = u./sqrt(4*u.^2 - 1);

    A = zeros(n);
    A(2:n+1:n*(n-1)) = u;
    A(n+1:n+1:n^2-1) = u;

    [v, x] = eig(A);
    [x, k] = sort(diag(x));    
    qpts = [x'; 2*v(1,k).^2];
end
