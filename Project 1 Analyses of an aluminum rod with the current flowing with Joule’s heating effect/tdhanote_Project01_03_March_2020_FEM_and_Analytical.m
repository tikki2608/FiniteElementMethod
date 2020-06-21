function [fe,K,T,I_max,e_L2] = tdhanote_Project01_03_March_2020_FEM_and_Analytical(nn_per_element,number_of_elements,k,r0, r1,L,rho,no_qpts)
%% input defination 

% nn_per_element = Number of nodes per element
% number_of_elements = Number of elements
% k = Thermal conductivity
% r0 = Radius at the ends of bar
% r1 = radius at the center of bar
% roh = Electrical Resistivity
% no_qpts = Number of quadrature Points

%% Inpts
if nargin==0
nn_per_element=3;  
number_of_elements=4;
k=205;
r0=0.002;
r1=0.001;
L=0.01;
I=1000;
rho=2.82e-8;
no_qpts=10;
end
count=1;
for number_of_elements=[4,8,16,32,64,128,256,512]
%% Heat Generation as a function of x

s=@(x)(((I^2)*rho)/(pi()*(r1 + (r0-r1)*(x/L))^2));
%% Arear as a function of x

A=@(x)(pi()*(r1 + (r0-r1)*(x/L))^2);

%% Quadture points selection

qpts=quadrature(no_qpts);

%% Total nodes and connectivity matrix determination Corresponding shape function
%Linear
if nn_per_element==2
    nodes_total=number_of_elements+1;
    conn=[1:number_of_elements;2:number_of_elements+1];
    shape=@shape2;
% Quadratic    
elseif nn_per_element==3
     nodes_total=2*number_of_elements+1;
     conn=[1:2:2*number_of_elements-1;2:2:2*number_of_elements;3:2:2*number_of_elements+1];
     shape=@shape3;
%Cubic     
else
 conn=[1:3:3*number_of_elements-2;2:3:3*number_of_elements-1;3:3:3*number_of_elements;4:3:3*number_of_elements+1];
 nodes_total=3*number_of_elements+1;
 shape=@shape4;
end
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
fe=zeros(nodes_total,1);

for c=conn
    xe=x(:,c);
        for i=1:no_qpts   
    [N,DNDp]=shape(qpts(1,i));
    J=xe*DNDp;
    w=qpts(2,i);
    xp=xe*N;
    F=s(xp)*N*J*w;
    fe(c)=fe(c)+F;
        end
        
end
%% Temperature Calculation at nodes

K(nodes_total, :)= 0;
K(nodes_total,nodes_total)=1;
fe(nodes_total,1)=20;
T=K\fe;

DTDX=[];
x1=[];
Heat_flux(1)=0;
%% Heat Flux 
for c=conn
    xe=x(:,c);
    for p=linspace(-1,1,nn_per_element) 
        [N,dNdp]=shape(p);
        J=xe*dNdp;
        B=dNdp/J;
        DTDX(end+1)=-B'*T(c);  
        x1(end+1)=xe*N;
    end    
end
Heat_flux(1:(size(DTDX,2)))=DTDX*k;


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
c2_I_max=solve(Ta(r0,r1,k,Ia,rho,L,x_T,c1_I_max,C2)==20,C2);
I_max=solve(Ta(r0, r1, k, Ia, rho, L, x_dTdx, c1_I_max, c2_I_max)==660,Ia);
I_max=max(double(I_max));
c1=solve(dTdxa(r0,r1,k,I,rho,L,x_dTdx,c1)==0,c1);
c2=solve(Ta(r0,r1,k,I,rho,L,x_T,c1,c2)==20,c2);
x_a=linspace(0,0.01,11);
Temp_analytical=Ta(r0,r1,k,I,rho,L,x_a,c1,c2);
Heat_flux_analytical=-HF(r0,r1,k,I,rho,L,x_a,c1);

% %% Plots
%  Abaqus_temp=xlsread('Abaqus_Linear_Element_Temperature.xlsx');
%  Abaqus_Heat_flux=xlsread('Abaqus_Linear Element_Heat_Flux.xlsx');
% figure('name',"Temprature Profile");
% plot(x,T,'o--')
% hold on
% plot(x_a,Temp_analytical)
% hold on
% plot(Abaqus_temp(1:end,1),Abaqus_temp(1:end,2),'b*--')
% xlim([0 L])
% title('Temperature Profile')
% xlabel('X Distance')
% ylabel('Temeperature')
% legend('FEM','Analytical','Abaqus')
% hold off
% figure('name',"Heat Flux Profile");
% % XX=linspace(0,L,length(Heat_flux));
% plot(x1,Heat_flux)
% hold on
% plot(x_a,Heat_flux_analytical)
% hold on
% plot(Abaqus_Heat_flux(1:end,1),Abaqus_Heat_flux(1:end,2),'b*--')
% xlim([0 L])
% title('Heat Flux Profile')
% xlabel('X distance')
% ylabel('Heat Flux')
% legend('FEM','Analytical','Abaqus')
% hold off

%% Convergence
Fun_num=0;
Fun_den=0;
for c=conn
    xe=x(c);
for q=1:size(qpts,2)
    [N,dNdp]=shape(qpts(1,q));
    xee=xe*N;
    T1=T(c)'*N;
    J=xe*dNdp;
    Temp_Exact=Ta(r0,r1,k,I,rho,L,xee,c1,c2);
    fun_num=(Temp_Exact-T1)^2;
    fun_den=Temp_Exact^2;
    Fun_num=Fun_num+fun_num*J*qpts(2,q);
    Fun_den=Fun_den+fun_den*J*qpts(2,q);
end
end
e_L2(count)=sqrt(abs(Fun_num/Fun_den));
LE(count)=L/number_of_elements;
Fun_s_num=0;
Fun_s_den=0;
 for c=conn
     xe=x(:,c);
     for p=linspace(-1,1,nn_per_element) 
         [N,dNdp]=shape(p);
         XX=xe*N;
         J=xe*dNdp;
         B=dNdp/J;  
         Heat_Exact=-HF(r0,r1,k,I,rho,L,XX,c1);
         heat_flux_FEM=-B'*T(c)*k;
         fun_s_num=(Heat_Exact-heat_flux_FEM)^2;
         fun_s_den=Heat_Exact^2;
         Fun_s_num=Fun_s_num+(fun_s_num*J*qpts(2,q))/2;
         Fun_s_den=Fun_s_den+(fun_s_den*J*qpts(2,q))/2;
     end
 end
  e_en(count)=sqrt(abs(Fun_s_num/Fun_s_den));
count=count+1;
end
figure('name',"L2 temperature error norm vs element size")
plot(log10(LE),log10(e_L2),'o')
title('log(eL2)vs log(LE)')
xlabel('log(LE)--log of element length')
ylabel('log(eL2)--log of L2 Temperature error norm')
legend('Quadratic Element')
figure('name',"flux error norm vs element size")
plot(log10(LE),log10(e_en),'o')
title('log(een)vs log(LE)')
xlabel('log(LE)--log of element length')
ylabel('log(een)--log of flux error norm')
legend('Quadratic Element')

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
