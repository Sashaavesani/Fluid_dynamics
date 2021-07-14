% es6
%
% Solve the diffeq problem
% { fu' + 2u" = 0
% { f' - u = 0
% { f(0) = 0; f'(inf) = 1
% { u(0) = 0; u(inf) = 1
% using discretisation with finite difference

clear all
close all


% problem's data
U = 3;                   % Exterior velocity
nu = 1.5;                % viscosity

% boundary conditions 
f1 = 0;
der_fm = 1;
u1 = 0;
um = 1;

% x and y range boundarys
xa = 0;
xb = 2; 
y1 = 0;
ym = 10;


m = 200;

% nodes 
y = zeros(m,1);
for i = 1:m
    y(i) = (ym-y1) * ((i-1)/(m-1))^2; 
end
h = diff(y);

eta = y * sqrt(U/(nu*xb));


d = 1./(h(1:m-2)+h(2:m-1));
F1 = spdiags([[-d;0;0],[0;0;d]],[-1,1],m,m);

U1 = spdiags([[-d;0;0],[0;0;d]],[-1,1],m,m);

d1 = 1./(h(1:m-2).*(h(1:m-2)+h(2:m-1)));
d2 = 1./(h(2:m-1).*(h(1:m-2)+h(2:m-1)));
U2 = spdiags([[2*d1;0;0],[0;-2*d1-2*d2;0],[0;0;2*d2]],[-1,0,1],m,m);
b = [zeros(2*m,1)];


% boundary conditions
U2(1,1) = 1/2;
U2(m,m) = 1/2;
b(m) = 1;
F1(m,m-1) = -1./h(end);     
F1(m,m) = 1./h(end);      
F1(1,1) = 1;


F = @(f,u) [f.*(U1*u) + 2*(U2*u) - b(1:m);...
    F1*f - u - b(m+1:2*m)];
JF = @(f,u) [diag(U1*u), f.*U1 + 2*U2;...
    F1, -eye(m)];

%_____________________________________ Newton
tol = 1e-8;
f0 = [0; ones(m-1,1)];  
u0 = [0; ones(m-1,1)];
sol0 = [f0;u0];
sol = sol0;
f = sol(1:m);
u = sol(m+1:2*m);
delta = -JF(f,u)\F(f,u);
 while (norm (delta, inf) > tol)
    sol = sol + delta;
    f = sol(1:m);
    u = sol(m+1:2*m);
    delta = -JF(f,u)\F(f,u);
 end
sol = sol+delta;
f = sol(1:m);
u = sol(m+1:2*m);

u_2 = U * u;
v_2 = U * 1/2 *sqrt(nu/U) * 1/sqrt(xb) * (eta.*u - f);

%_______________________________________ plot 
figure
plot(u_2, y, v_2, y, 'LineWidth',1.5);
axis([0 3 0 10]);


%% es5

% boundary conditions 
u_x0 = 0;
v_x0 = 0;
u_x10 = U;

% number of x e y nodes 
mx = 50;
my = 100;

% initial conditions
u_0y = U.*ones(my,1);
v_0y = zeros(my,1);


%% x nodes
x = zeros(mx,1);
for i = 1:mx
    x(i) = (xb-xa) * ((i-1)/(mx-1))^2; 
end
hx = diff(x);

%% y nodes 
y = zeros(my,1);
for i = 1:my
    y(i) = (ym-y1) * ((i-1)/(my-1))^2; 
end
hy = diff(y);

D1v = spdiags([[0;1./hy(1:end);0], [-1./hy(1:end);0;0]], [0,-1], my,my);

d1 = 1./(hy(1:my-2).*(hy(1:my-2)+hy(2:my-1)));
d2 = 1./(hy(2:my-1).*(hy(1:my-2)+hy(2:my-1)));
D2 = spdiags([[2*d1;0;0],[0;-2*d1-2*d2;0],[0;0;2*d2]],[-1,0,1],my,my);

d = 1./(hy(1:my-2)+hy(2:my-1));
D1uv = spdiags([[-d;0;0],[0;0;d]],[-1,1],my,my);

D1x_i = spdiags([[1/2*ones(my-1,1);0;0], [0;1/2*ones(my-1,1);0]],[-1,0],my,my);
D1x_iprev = D1x_i;

D1x_isquare = speye(my);
D1x_isquareprev = D1x_isquare;

D1v(1,1) = 1;
D1x_i(end,:) = zeros(1,my);
D1x_iprev(end,:) = zeros(1,my);
D2(1,1) = 1;
D2(end,end) = 1;
D1x_isquare(1,:) = zeros(1,my);
D1x_isquareprev(1,:) = zeros(1,my);
D1x_isquare(end,:) = zeros(1,my);
D1x_isquareprev(end,:) = zeros(1,my);

b = [zeros(my-1,1);-nu*U];


u_Eprevious = [0;u_0y(2:my)];                   % first previous point for Euler 
                                                % (first component =0 to sadisfy bound. cond.)
v = v_0y;                                       % initial condition of v

for i = 1:mx-1
    
    F = @(u,v) [1/hx(i).*(D1x_i*u) - 1/hx(i).*(D1x_iprev*u_Eprevious) + D1v*v; ...
        1/hx(i).*(D1x_isquare*(u.^2)) - 1/hx(i).*(D1x_isquareprev*(u_Eprevious.^2)) + D1uv*(u.*v) - nu*(D2*u) - b];
    JF = @(u,v) [1/hx(i).*(D1x_i), D1v; ...
        2/hx(i).*diag(D1x_isquare*u) + D1uv*diag(v) - nu*D2, D1uv*diag(u)];
    
    %____________________________________ Newton
    u0 = u_Eprevious;
    v0 = v;
    % let's initial points satisfy boundary conditions
    u0(1) = 0;
    v0(1) = 0;
    u0(end) = U;
    sol0 = [u0;v0];
    sol = sol0;
    u = sol(1:my);
    v = sol(my+1:2*my);
    
    delta = -JF(u,v)\F(u,v);
    tol = 1e-6;
    maxits = 20;
    its = 1;
     while (norm (delta, inf) > tol) && (its<maxits)
        sol = sol + delta;
        u = sol(1:my);
        v = sol(my+1:2*my);
        delta = -JF(u,v)\F(u,v);
        its = its + 1;
     end
    sol = sol+delta;
    u = sol(1:my);                  % it will became the next u_Eprevious
    v = sol(my+1:2*my);
    u_Eprevious = u;   
end 
hold on
plot(u,y,'o',v,y,'o')
legend('ues6','ves6','ues5','ves5')
