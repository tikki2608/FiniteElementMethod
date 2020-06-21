function [Temperature,HeatFlux]= tdhanote_project_Manufactured_Solution_29_April_2020_Revision_2(h,r0,r1,etype,R,sbc,nqpts,T0,k)
if nargin == 0
    h=1;
    r0=10;
    r1=20;
    etype='t6';     % q4 q8 q9 t3 t6
    nqpts=5;
    k=0.017;           %W/mm-K
    th=1;
end
%% Manufactured Solution Method
H=[2*h h h/2 h/4 h/8];
T_star=@(X,Y) (sin(X) + sin(Y));
Q_heat=@(X,Y) ([-k*cos(X) , -k*cos(Y)]);
%% Heat source %%
S_heat=@(X,Y) (k*(sin(X) + sin(Y)));
e_L2=[];
e_en=[];
leng=[];
for h=H
%% Shape function and quarature points
[mesh] = make_project2_mesh(h, r0, r1, etype);
[Shape,qpts,facep,shape_edge,nodes_outer,edgeconn_inner]=SQF(etype,nqpts,mesh,r0,r1);
%% mesh %%
x=mesh.x;
conn=mesh.conn;
%% Average Temperture calculation
K=spalloc(length(x),length(x),2*length(x));
Fe=zeros(length(x),1);
    for c=conn
        xe=mesh.x(:,c);
        Ke=zeros(length(c));
        for q=qpts
            [N,dNdp]=Shape(q);
            J=xe*dNdp;
            B=dNdp/J;
            w=q(end);
            Ke=Ke+B*k*eye(2)*th*B'*det(J)*w;
            f=S_heat(xe(1,:)*N,xe(2,:)*N);
            Fe(c)=Fe(c)+N*f*det(J)*w;
        end
        K(c,c)=K(c,c)+Ke;
    end
    for ec=edgeconn_inner
        xe=mesh.x(:,ec);
        for q=quadrature_1D(nqpts)
            [N,dNdp]=shape_edge(q(1));
            J=xe*dNdp;
            n=[0 -1;1 0]*J/norm(J);
            XE=xe*N;
            Q=dot(n,Q_heat(XE(1),XE(2)));
            Fe(ec)=Fe(ec)-N*Q*norm(J)*q(end);      
        end
    end
    fixed_nodes_temp=nodes_outer(1:end);
    K(fixed_nodes_temp, :)= 0;
    K(fixed_nodes_temp,fixed_nodes_temp)=eye(length(fixed_nodes_temp));
    XY=x(:,fixed_nodes_temp);
    Fe(fixed_nodes_temp,1)=(T_star(XY(1,:),XY(2,:)))';
    T=K\Fe;
    %% Convergence
Fun_num=0;
Fun_den=0;
for c=conn
    xe=x(:,c);
    len=norm(xe(:,2)-xe(:,1));
for q=qpts
    [N,dNdp]=Shape(q(1:2));
    xee=xe*N;
    T1=T(c)'*N;
    J=xe*dNdp;
    Temp_Exact=(T_star(xee(1),xee(2)));
    fun_num=(Temp_Exact-T1)^2;
    fun_den=Temp_Exact^2;
    Fun_num=Fun_num+fun_num*det(J)*q(end);
    Fun_den=Fun_den+fun_den*det(J)*q(end);
end
end

Fun_s_num=0;
Fun_s_den=0;
for c=conn
    xe=x(:,c);
    for q=qpts
        [N,dNdp]=Shape(q(1:2));
        J=xe*dNdp;
        B=dNdp/J;
        Tf=T(c);
        Xe=xe*N;
        n=Xe/norm(Xe);
        Q=-dot(n,Tf'*B*eye(2)*k);
        Heat_Exact=dot(n,Q_heat(Xe(1),Xe(2)));
        fun_s_num=(Heat_Exact-Q)^2;
        fun_s_den=Heat_Exact^2;
        Fun_s_num=Fun_s_num+(fun_s_num*det(J)*q(end))/2;
        Fun_s_den=Fun_s_den+(fun_s_den*det(J)*q(end))/2;
    end
end
e_L2(end+1)=sqrt(abs(Fun_num/Fun_den));
e_en(end+1)=sqrt(abs(Fun_s_num/Fun_s_den));
leng(end+1)=len;
end
figure('name',"log(h) vs log(e_L2)");
plot(log(H),log(e_L2),'--o')
title(['log(eL2) vs log(h) for ',etype ,' Element'])
xlabel('log(h)--log of element length')
ylabel('log(eL2)--log of L2 Temperature error norm')
legend([etype ,' Element'])
figure('name',"log(h) vs log(e_en)");
plot(log(H),log(e_en),'--o')
title(['log(een)vs log(LE)for ',etype ,' Element'])
xlabel('log(h)--log of element length')
ylabel('log(een)--log of flux error norm')
legend([etype ,' Element'])
end


function[N, dNdp] = shape_q4(p)
N = [(1/4)*(1-p(1))*(1-p(2));
    (1/4)*(1+p(1))*(1-p(2));
    (1/4)*(1+p(1))*(1+p(2));
    (1/4)*(1-p(1))*(1+p(2))];
dNdp = [(1/4)*(p(2)-1),(1/4)*(p(1)-1);
    (1/4)*(1-p(2)),(-1/4)*(p(1)+1);
    (1/4)*(1+p(2)),(1/4)*(1+p(1));
    (-1/4)*(1+p(2)),(1/4)*(1-p(1))];
end
function[N, dNdp] = shape_q8(p) 
N = [(-1/4)*(1-p(1))*(1-p(2))*(1+p(1)+p(2));    %1-1
    (-1/4)*(1+p(1))*(1-p(2))*(1-p(1)+p(2));     %3-2
    (-1/4)*(1+p(1))*(1+p(2))*(1-p(1)-p(2));     %5-3
    (-1/4)*(1-p(1))*(1+p(2))*(1+p(1)-p(2));     %7-4
    (1/2)*(1-p(1))*(1+p(1))*(1-p(2));           %2-5
    (1/2)*(1+p(1))*(1+p(2))*(1-p(2));           %4-6
    (1/2)*(1-p(1))*(1+p(1))*(1+p(2));           %6-7
    (1/2)*(1-p(1))*(1+p(2))*(1-p(2))];          %8-8
dNdp = [(-1/4)*(p(2)-1)*(2*p(1)+p(2)),(-1/4)*(p(1)-1)*(2*p(2)+p(1));    %1-1
     (1/4)*(p(2)-1)*(-2*p(1)+p(2)),(1/4)*(p(1)+1)*(2*p(2)-p(1));        %3-2
     (1/4)*(1+p(2))*(2*p(1)+p(2)),(1/4)*(1+p(1))*(p(1)+2*p(2));         %5-3
     (-1/4)*(1+p(2))*(p(2)-2*p(1)),(-1/4)*(p(1)-1)*(2*p(2)-p(1));       %7-4
     p(1)*(p(2)-1),(1/2)*(1+p(1))*(-1+p(1));                            %2-5 
     (-1/2)*(p(2)+1)*(p(2)-1),(-1)*(1+p(1))*p(2);                       %4-6
     (-1)*p(1)*(1+p(2)),(-1/2)*(1+p(1))*(p(1)-1);                       %6-7
     (1/2)*(1+p(2))*(p(2)-1),p(2)*(p(1)-1)];                            %8-8
end
function[N, dNdp] = shape_q9(p)
N = [0.5*p(1)*(p(1)-1)*0.5*p(2)*(p(2)-1);
    0.5*p(1)*(p(1)+1)*0.5*p(2)*(p(2)-1);
    0.5*p(1)*(p(1)+1)*0.5*p(2)*(p(2)+1);
    0.5*p(1)*(p(1)-1)*0.5*p(2)*(p(2)+1);
    (1-p(1)^2)*0.5*p(2)*(p(2)-1);
    0.5*p(1)*(p(1)+1)*(1-p(2)^2);
    (1-p(1)^2)*0.5*p(2)*(p(2)+1);
    0.5*p(1)*(p(1)-1)*(1-p(2)^2);
    (1-p(1)^2)*(1-p(2)^2)];
dNdp = [0.5*(2*p(1)-1)*0.5*p(2)*(p(2)-1),0.5*p(1)*(p(1)-1)*0.5*(2*p(2)-1);
    0.5*(2*p(1)+1)*0.5*p(2)*(p(2)-1),0.5*p(1)*(p(1)+1)*0.5*(2*p(2)-1);
    0.5*(2*p(1)+1)*0.5*p(2)*(p(2)+1),0.5*p(1)*(p(1)+1)*0.5*(2*p(2)+1);
    0.5*(2*p(1)-1)*0.5*p(2)*(p(2)+1),0.5*p(1)*(p(1)-1)*0.5*(2*p(2)+1);
    (-2*p(1))*0.5*p(2)*(p(2)-1),(1-p(1)^2)*0.5*(2*p(2)-1);
    0.5*(2*p(1)+1)*(1-p(2)^2),0.5*p(1)*(p(1)+1)*(-2*p(2));
    (-2*p(1))*0.5*p(2)*(p(2)+1),(1-p(1)^2)*0.5*(2*p(2)+1);
    0.5*(2*p(1)-1)*(1-p(2)^2),0.5*p(1)*(p(1)-1)*(-2*p(2));
    (-2*p(1))*(1-p(2)^2),(1-p(1)^2)*(-2*p(2))];
end

function[N, dNdp] = shape_t3(p)
N = [p(1);
    p(2);
    1 - p(1) - p(2)];
dNdp = [1, 0;
    0, 1;
    -1, -1];
end
function[N, dNdp] = shape_t6(p)
N = [p(1)*(2*p(1)-1);
    p(2)*(2*p(2)-1);
    (1-p(2)-p(1))*(2*(1-p(2)-p(1))-1);
    4*p(1)*p(2);
    4*p(2)*(1-p(2)-p(1));
    4*p(1)*(1-p(2)-p(1))];
    
dNdp =  [4*p(1) - 1, 0 ;
        0, 4*p(2) - 1;
        -3+4*p(1)+4*p(2),-3+4*p(1)+4*p(2);
        4*p(2), 4*p(1);
        -4*p(2),(4-8*p(2)-4*p(1));
        4-4*p(2)-8*p(1),-4*p(1)];       
end
function [N,dNdp] = shape2(p)
N=0.5*[1-p(1);1+p(1)];
dNdp=[-0.5;0.5];
end
function [N, dNdp] = shape3(p) % Define quadratic shape functions and derivatives 
N=1/2*[-(1-p)*p ; 2*(1-p)*(1+p); (1+p)*p ];
dNdp=[ p-0.5000 ; -2*p ; p+0.5000];
end

%% Gaussian Quadrature 1D 
function [qpts] = quadrature_1D(n)
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
function [Shape,qpts,facep,shape_edge,nodes_outer,edgeconn_inner]=SQF(etype,nqpts,mesh,r0,r1)
x=mesh.x;
edge_inner=find(abs(x(1,:).^2+x(2,:).^2-r0^2)<1e-10);
edge_outer=find(abs(x(1,:).^2+x(2,:).^2-r1^2)<1e-10);
edge_sidex=find(abs(x(2,:))<1e-10);
edge_sidey=find(abs(x(1,:))<1e-10);
nodes_outer=[edge_sidex edge_outer(2:end-1) edge_sidey];
if strcmp(etype, 'q4')
    Shape=@shape_q4;
    shape_edge=@shape2;
    qpts=Quadrature_2D_Quadilateral_element(nqpts);
    facep=mesh.conn';
    edgeconn_inner=[edge_inner(1:end-1);edge_inner(2:end)];
elseif strcmp(etype, 't3')
    Shape=@shape_t3;
    shape_edge=@shape2;
    qpts=   [0.1012865073 0.1012865073 0.0629695903;
        0.7974269853 0.1012865073 0.0629695903;
        0.1012865073 0.7974269853 0.0629695903;
        0.4701420641 0.0597158717 0.0661970764;
        0.4701420641 0.4701420641 0.0661970764;
        0.0597158717 0.4701420641 0.0661970764;
        0.3333333333 0.3333333333 0.1125]';
    facep=mesh.conn';
    edgeconn_inner=[edge_inner(1:end-1);edge_inner(2:end)];
elseif strcmp(etype, 'q8')
    Shape=@shape_q8;
     shape_edge=@shape3;
    qpts=Quadrature_2D_Quadilateral_element(nqpts);
    facep=mesh.pconn';
    edgeconn_inner=[edge_inner(1:2:end-2);edge_inner(2:2:end-1);edge_inner(3:2:end)];
elseif strcmp(etype, 'q9')
    Shape=@shape_q9;
    qpts=Quadrature_2D_Quadilateral_element(nqpts);
    facep=mesh.pconn';
elseif strcmp(etype, 't6')
    Shape=@shape_t6;
     shape_edge=@shape3;
    qpts=   [0.1012865073 0.1012865073 0.0629695903;
        0.7974269853 0.1012865073 0.0629695903;
        0.1012865073 0.7974269853 0.0629695903;
        0.4701420641 0.0597158717 0.0661970764;
        0.4701420641 0.4701420641 0.0661970764;
        0.0597158717 0.4701420641 0.0661970764;
        0.3333333333 0.3333333333 0.1125]';
    facep=mesh.pconn';
    edgeconn_inner=[edge_inner(1:2:end-2);edge_inner(2:2:end-1);edge_inner(3:2:end)];
end
end