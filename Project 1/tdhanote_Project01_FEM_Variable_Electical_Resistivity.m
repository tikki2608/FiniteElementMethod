function [fe,K,T,I_max] = tdhanote_Project01_FEM_Variable_Electical_Resistivity(nn_per_element,number_of_elements,k,r0, r1,L,T0,no_qpts)
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
nn_per_element=4;  
number_of_elements=1;
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
dTT=[];
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
%% convergence of Temperature for Electrical resistivity change
dTT(:,count+1)=(dT./T);
count=count+1;
end

%% Plots
figure('name',"Temprature Profile");
plot(x,T,'o--')
xlim([0 L])
title('Temperature Profile')
xlabel('X Distance (m)')
ylabel('Temeperature(^0C)')
legend('FEM')
hold off
figure('name',"Convergence plot Temperature each Node");
plot(1:count,dTT,'o--')
title('Convergence plot Temperature each Node')
xlabel('Number of Iteration')
ylabel('Error')
legend('Node 1','Node 2','Node 3','Node 4')
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
