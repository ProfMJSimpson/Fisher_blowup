%Implicit NR solver for Fisher-Stefan on a geometrically-spaced mesh
clearvars
N=2000; %N is the number of intervals
dxi0 = 1e-6; %smallest grid space on RHS
ff = @(x) (1-x^N)/(1-x)-1/dxi0;
ee = fsolve(ff,1.01);
% format long
% ee
xi=zeros(1,N+1);
dxi=zeros(1,N);
xi(1,N+1)=1.0;
term=0;
for i=N:-1:1
    term=term+1;
    xi(1,i) = xi(1,i+1)-dxi0*ee^(term-1);
    dxi(1,i) = xi(1,i+1)-xi(1,i);
end
xi(1,1)=0;
figure 
plot(dxi,'.')
 xlabel('N')
 ylabel('dxi')
 
 figure
 plot(dxi,'.')
 xlabel('N')
 ylabel('dxi')
 set(gca, 'YScale', 'log')

dt=0.0001;
T=6.5;
tol=1e-10;
L0=1000.0;
l=L0;
pl=l;
kappa=-1.01;
maxsteps=round(T/dt);
Lrecord=zeros(maxsteps+1,1);
Lrecord(1,1)=L0;
wsrecord=zeros(maxsteps+1,1);
wsrecord(1,1)=0;



a=zeros(1,N+1);
b=zeros(1,N+1);
c=zeros(1,N+1);
d=zeros(1,N+1);
u=0.5*ones(1,N+1);
pu=0.5*ones(1,N+1);
delu=ones(1,N+1);
u(1,N+1)=0;
pu(1,N+1)=0;


t=0.0;


for i=1:maxsteps
    t=t+dt;
    kk=0;
    delu=ones(1,N+1);
    while norm(delu,Inf) > tol
        kk=kk+1;
        
a(1,1)=0.0;
b(1,1)=-1.0;
c(1,1)=1.0;
d(1,1)=-1*(u(1,2)-u(1,1));

a(1,N+1)=0.0;
b(1,N+1)=1.0;
c(1,N+1)=0.0;
d(1,N+1)=0.0-u(1,N+1);

for j=2:N
    hp=dxi(j);
    hm=dxi(j-1);
    alpha = 1/(hm*(hm+hp));
    beta = -1/(hm*hp);
    delta = 1/(hp*(hm+hp));
    
a(1,j)=2*alpha/l^2-xi(1,j)*(l-pl)*alpha*hp/(l*dt);
b(1,j)=-1.0/dt+1.0-2.0*u(1,j)+2*beta/l^2-xi(1,j)*(l-pl)*beta*(hm-hp)/(l*dt);
c(1,j)=2*delta/l^2+xi(1,j)*(l-pl)*delta*hm/(l*dt);
d(1,j)=(u(1,j)-pu(1,j))/dt-2*(alpha*u(1,j-1)+beta*u(1,j)+delta*u(1,j+1))/l^2 ...
-xi(1,j)*(l-pl)*(delta*hm*u(1,j+1)+beta*(hm-hp)*u(1,j)-alpha*hp*u(1,j-1))/(l*dt)...
-u(1,j)*(1-u(1,j));    

end

delu = thomas(N+1,a,b,c,d);

u(1,:)=u(1,:)+delu(1,:); 

%l=pl+dt*kappa*u(1,N)/(pl*dxi0);

    hp=dxi(N);
    hm=dxi(N-1);
    alpha = 1/(hm*(hm+hp));
    beta = -1/(hm*hp);

l=pl-(dt*kappa/(pl))*(alpha*u(1,N-1)*hp+beta*u(1,N)*(hp+hm));
    end
    
    if mod(i,1000)==0
    fprintf('Time %d\n',t); 
    fprintf('Iteration %d\n',kk);
    end
    
Lrecord(i+1,1)=l; 
wsrecord(i+1,1)=(l-pl)/dt;
pu(1,:)=u(1,:);
pl=l;
Lrecord(i+1,1)=l;
end

for i=1:N+1
    xx(1,i)=l*xi(1,i);
end
figure
plot(xx,u)
axis([0 1.1*l 0 1.1])
figure
plot(dt*(1:maxsteps+1),Lrecord)


figure
plot(dt*(1:maxsteps+1),wsrecord)
maxspeed = max(abs(wsrecord));
%ylim([-maxspeed,maxspeed]);
xlim([0.1,30]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Algorithm
function x = thomas(N,a,b,c,d)
x=zeros(1,N);
    bb=b;
    dd=d;
    for i=2:N
        ff=a(i)/bb(i-1);
        bb(i)=bb(i)-c(i-1)*ff;
        dd(i)=dd(i)-dd(i-1)*ff;
    end
    
    for i=1:N-1
    x(N)=dd(N)/bb(N);    
    j=N-i;
    x(j)=(dd(j)-c(j)*x(j+1))/bb(j);
    end
end



