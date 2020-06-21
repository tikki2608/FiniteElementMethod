function [T, Heat_flux] = tdhanote_Project01_FEM(nn_per_element,number_of_elements,k,r0, r1,L,rho,no_qpts)
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
no_qpts=20;
end
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
%% Plots
figure('name',"Temprature Profile");
plot(x,T)
hold on
xlim([0 L])
title(['Temperature Profile with ',num2str(number_of_elements), ' Cubic elements'])
xlabel('x Distance (m)')
ylabel('Temeperature (^0C)')
legend([num2str(number_of_elements), ' Cubic elements'])
hold off
figure('name',"Heat Flux Profile");
plot(x1,Heat_flux)
xlim([0 L])
title(['Heat Flux Profile with ',num2str(number_of_elements), ' Cubic elements'])
xlabel('x Distance (m)')
ylabel('Heat Flux (W/m^2)')
legend([num2str(number_of_elements), ' Cubic elements'])
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