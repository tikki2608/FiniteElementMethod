function [Temperature,HeatFlux]= tdhanote_project_02_27_April_2020_Revision_8(h,r0,r1,etype,R,sbc,nqpts,T0,k)
if nargin == 0
    h=1;
    r0=10;          %milimeter
    r1=20;          %milimeter
    etype='q4';     % q4 q8 q9 t3 t6
    R=3;            %milimeter
    sbc=5.670e-14;   %Wmm-2K-4
    nqpts=5;
    T0=300;         % K
    k=eye(2)*0.017;           %W/mm-K
    th=1;
end
%% input variables %%
xc=0.5*(r1+r0)*cos(pi/4);
yc=0.5*(r1+r0)*sin(pi/4);
Ac=pi*(r1^2-r0^2)/4;
%% Shape function and quarature points
[mesh] = make_project2_mesh(h, r0, r1, etype);
if strcmp(etype, 'q4')
    Shape=@shape_q4;
    qpts=Quadrature_2D_Quadilateral_element(nqpts);
    facep=mesh.conn';
elseif strcmp(etype, 't3')
    Shape=@shape_t3;
    qpts=   [0.1012865073 0.1012865073 0.0629695903;
        0.7974269853 0.1012865073 0.0629695903;
        0.1012865073 0.7974269853 0.0629695903;
        0.4701420641 0.0597158717 0.0661970764;
        0.4701420641 0.4701420641 0.0661970764;
        0.0597158717 0.4701420641 0.0661970764;
        0.3333333333 0.3333333333 0.1125]';
    facep=mesh.conn';
elseif strcmp(etype, 'q8')
    Shape=@shape_q8;
    qpts=Quadrature_2D_Quadilateral_element(nqpts);
    facep=mesh.pconn';
elseif strcmp(etype, 'q9')
    Shape=@shape_q9;
    qpts=Quadrature_2D_Quadilateral_element(nqpts);
    facep=mesh.pconn';
elseif strcmp(etype, 't6')
    Shape=@shape_t6;
    qpts=   [0.1012865073 0.1012865073 0.0629695903;
        0.7974269853 0.1012865073 0.0629695903;
        0.1012865073 0.7974269853 0.0629695903;
        0.4701420641 0.0597158717 0.0661970764;
        0.4701420641 0.4701420641 0.0661970764;
        0.0597158717 0.4701420641 0.0661970764;
        0.3333333333 0.3333333333 0.1125]';
    facep=mesh.pconn';
end
%% Heat source %%
s=@(x,y) (exp(-(((x-xc)^2+(y-yc)^2)/R^2)));
%% mesh %%
x=mesh.x;
conn=mesh.conn;
%% Average Temperture calculation
fe=0;
for c=conn
    xe=x(:,c);
    for q=qpts
        [N,dNdp]=Shape(q(1:2));
        J=xe*dNdp;
        xp=xe*N;
        S=s(xp(1),xp(2));
        F=S*det(J)*q(end);
        fe=fe+F;
    end
end
Tavg=((fe/(Ac*sbc))+T0^4)^(1/4);
Heff=@(T) (sbc*(T^2+T0^2)*(T+T0));
K=zeros(length(x),length(x));
H=zeros(length(x),length(x));
T=zeros(length(x),1);
Ta=ones(length(x),1)*Tavg;
Fh=zeros(length(x),1);
Fe=zeros(length(x),1);
dT=100;
count=0;
while max(dT)>=0.01
    for c=conn
        xe=mesh.x(:,c);
        Ke=zeros(length(c));
        Hk=zeros(length(c));
        TA=Ta(c);
        for q=qpts
            [N,dNdp]=Shape(q);
            J=xe*dNdp;
            B=dNdp/J;
            w=q(end);
            Ke=Ke+B*k*th*B'*det(J)*w;
            f=s(xe(1,:)*N,xe(2,:)*N);
            TAP=TA'*N;
            hk=N*Heff(TAP)*N'*det(J)*w;
            Hk=Hk+hk;
            fh=N*Heff(TAP)*det(J)*w*T0;
            Fh(c)=Fh(c)+fh;
            Fe(c)=Fe(c)+N*f*det(J)*w;
        end
        K(c,c)=K(c,c)+Ke;
        H(c,c)=H(c,c)+Hk;
    end
    T=((K+H)\(Fe+Fh));
    dT=abs(T-Ta);
    Ta=T;
    count=count+1;
    plot(count,max(dT),'--ob')
    hold on
end
hold off
figure()
patch('Faces',facep,'Vertices',mesh.x','FaceVertexCData',T,'FaceColor','interp');
colorbar
%%  conn and pconn need to be used in this case.. (works now for both temperature and flux plots if u want to use)
%% better looking than above code as lines are thinner 
% p.vertices =mesh.x';
% p.faces = mesh.conn';
% p.facecolor ='interp';
% p.facevertexcdata=T; %plot temp
% p.edgealpha = 0.2; %transparency
% clf;
% patch(p);
% colorbar
%% flux
%
%% [Add this to the above code once we get correct answer] %%
%
SpA = spalloc(length(mesh.x), length(mesh.x), 9*(length(mesh.x)));
den = zeros(length(mesh.x),2);
for c=conn
    xe=x(:,c);
    SpAe = zeros(length(c));
    for q=qpts
        [N,dNdp]=Shape(q);
        J=xe*dNdp;
        B=dNdp/J;
        w=q(end);
        Tf=T(c);
        Q=-Tf'*B*k;
        SpAe = SpAe+N*N'*det(J)*w;
        den(c,:) = den(c,:) + N*Q*det(J)*w;
    end
    SpA(c,c) = SpA(c,c) + SpAe;
end
flux=SpA\den;
Flux = sqrt(flux(:,1).^2+flux(:,2).^2);
figure()
patch('Faces',facep,'Vertices',mesh.x','FaceVertexCData',Flux,'FaceColor','interp');
colorbar

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