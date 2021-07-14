% ESERCIZIO 5 

% NUMERICAL SOLUTION OF BOUNDARY LAYER OF A BOUNDING SURFACE (2D)
% Solve the PDE problem:
% { ∂u/∂x + ∂v/∂y = 0;                      CONTINUITY EQUATION
% { ∂u^2/∂x + ∂uv/∂y = nu*∂u^2/∂y^2         MOMENTUM EQUATION
% { u(x,0) = v(x,0) = 0; u(x,inf) = U;      Boundary Conditions
% { u(0,y) = U; v(0,y) = 0;                 Initial Conditions

clear all 
close all 

tic
% problem's data
U = 3;                   % Exterior velocity
nu = 1.5;                % viscosity

% boundary conditions 
u_x0 = 0;
v_x0 = 0;
u_x10 = U;
dv_x10 = 0;

% number of x e y nodes 
mx = 100;
my = 200;

% initial conditions
u_0y = U.*ones(my,1);
v_0y = zeros(my,1);

% x and y range boundarys
xa = 0;
xb = 2; 
y1 = 0;
ym = 10;


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


%% FINITE DIFFERENCE IN y (RELATED ONLY TO INDEX j)
% We now create the matrices of finite difference for the terms divided by
% hy

%________________________ Matrix of finite difference for the term ∂v/∂y
% We use the following discretisation: (v_j-v_j-1)/(y_j-y_j-1)
% FIRST ROW OF D1v IS ZEROS

D1v = spdiags([[0;1./hy(1:end);0], [-1./hy(1:end);0;0]], [0,-1], my,my);


%________________________ Matrix of finite difference for the term ∂^2u/∂y^2
% We use the following discretization: [(u_j+1 - u_j)/(y_j+1 - y_j) - (u_j
% - u_j-1)/(y_j - y_j-1)] /[(y_j+1 - y_j-1)/2]
% FIRST ROW AND LAST ROW OF D2 ARE ZEROS

% denominator of subdiagonal terms
d1 = 1./(hy(1:my-2).*(hy(1:my-2)+hy(2:my-1)));
% denominator of superdiagonal terms 
d2 = 1./(hy(2:my-1).*(hy(1:my-2)+hy(2:my-1)));
D2 = spdiags([[2*d1;0;0],[0;-2*d1-2*d2;0],[0;0;2*d2]],[-1,0,1],my,my);


%________________________ Matrix of finite difference for the term ∂uv/∂y
% We use the following discretisation: (u_j+1*v_j+1 - u_j-1*v_j-1)/(y_j+1-y_j-1)
% FIRST ROW AND LAST ROW OF D1uv ARE ZEROS

d = 1./(hy(1:my-2)+hy(2:my-1));
D1uv = spdiags([[-d;0;0],[0;0;d]],[-1,1],my,my);


%% MATRICES FOR THE TERMS ∂u/∂x AND ∂u^2/∂x (RELATED ONLY TO INDEX j)
% REMARK: THESE DISCRETIZATIONS COME FROM EULER METHOD, i.e. ∂u/∂x = 
% (u_next-u_previous)/(x_next-x_previous)

%________________________ Matrices for the term ∂u/∂x (WE CONSIDER ONLY j-TERMS, so we 
% have (u_j+u_j-1)/2)
% FOR u_i
D1x_i = spdiags([[1/2*ones(my-1,1);0;0], [0;1/2*ones(my-1,1);0]],[-1,0],my,my);
% FOR u_i-1
D1x_iprev = D1x_i;

%________________________ Matrices for the term ∂u^2/∂x (WE CONSIDER ONLY j-TERMS, so we 
% have u_j.^2
% FOR u_i
D1x_isquare = speye(my);
% FOR u_i-1
D1x_isquareprev = D1x_isquare;


%% BOUNDARY CONDITIONS 
% from the fist equation of the system, through D1v we obtain the conditions
% for v
% from the second equation of the system, through D2 we obtain the conditions
% for u
% -----> AIM: the discretized system is:
% from line 1 to my
% { 1/hx(i).*(D1x_i*u) - 1/hx(i).*(D1x_iprev*u_Eprevious) + D1v*v = 0
% from line my+1 to 2my
% { 1/hx(i).*(D1x_isquare*(u.^2)) - 1/hx(i).*(D1x_isquareprev*(u_Eprevious.^2)) 
% + D1uv*(u.*v) - nu*(D2*u) = 0

% We want to obtain conditions for v with D1v
% _______LINE 1: considering that the first row of D1x_i and of D1x_iprev are zeros,
% we obtain D1v*v1 = 0. Let D1v(1,1) be 1, then v1 = 0 (FIRST BOUND. CONDITION)
% _______LINE my: changing the last row of D1x_i and D1x_iprev in zeros 
% ((D1x_i(end,:) = zeros(1,my); D1x_iprev(end,:) = zeros(1,my);)
% we get: D1v*v_my = 0. Remeber that D1v is the matrix of first derivative 
% of v so this equation tells us that the first derivative of v_my
% is 0 (and this is correct: far from the surface there 
% is only horizontal constant velocity U)
D1v(1,1) = 1;
D1x_i(end,:) = zeros(1,my);
D1x_iprev(end,:) = zeros(1,my);

% We want to obtain conditions for u with D2
% _______LINE my+1: changing the first row of D1x_isquare and D1x_isquareprev in zeros 
% (D1x_isquare(1,:) = zeros(1,my); D1x_isquareprev(1,:) = zeros(1,my);)
% we get: { D1uv*(u1.*v1) - nu*(D2*u1) = 0. Remember that the first row of
% D1uv is zeros so: {- nu*(D2*u1) = 0. It is necessary now to put D2(1,1) =
% 1 ---> {- nu*u1 = 0 ---> {u1 = 0 (SECOND BOUND. CONDITION)
% _______LINE 2*my+1: with the same argument we obtain {- nu*u_ym = 0. Now
% we want that u_ym = U, so let introduce b (known term) with b(end) =
% -nu*U. In details we have {- nu*u_ym = b(end) ---> {- nu*u_ym = -nu*U  
% ---> {u_ym = U (THIRD BOUND. CONDITION)

D2(1,1) = 1;
D2(end,end) = 1;

D1x_isquare(1,:) = zeros(1,my);
D1x_isquareprev(1,:) = zeros(1,my);
D1x_isquare(end,:) = zeros(1,my);
D1x_isquareprev(end,:) = zeros(1,my);

b = [zeros(my-1,1);-nu*U];


%% INITIAL CONDITION GIVEN BY PROBLEM
% REMEMBER: WE USE EULER METHOD in ∂u/∂x = -∂v/∂y. Infact for Euler we have
% ∂u/∂x = (u_i - u_i-1)/(x_i-x_i-1) so u_i = u_i-1 + (x_i-x_i-1)*∂u/∂x 
% = u_i-1 + (x_i-x_i-1)*(-∂v/∂y). 
% REMARK: IN CONTINUITY EQUATION WE USE INTERMEDIATED POINTS IN y.
% THEN FINITE DIFFERENCE METHOD WORKS ON BOTH u AND v.


u_Eprevious = [0;u_0y(2:my)];                   % first previous point for Euler 
                                                % (first component = 0 to satisfy bound. cond.)
v = v_0y;                                       % initial condition of v

for i = 1:mx-1
    % uE previous = aggiorno con eulero
    F = @(u,v) [1/hx(i).*(D1x_i*u) - 1/hx(i).*(D1x_iprev*u_Eprevious) + D1v*v; ...
        1/hx(i).*(D1x_isquare*(u.^2)) - 1/hx(i).*(D1x_isquareprev*(u_Eprevious.^2)) + D1uv*(u.*v) - nu*(D2*u) - b];
    JF = @(u,v) [1/hx(i).*(D1x_i), D1v; ...
        2/hx(i).*diag(D1x_isquare*u) + D1uv*diag(v) - nu*D2, D1uv*diag(u)];
    
    %____________________________________ Newton
    u0 = u_Eprevious;
    v0 = v;
    % let's initial points always satisfy boundary conditions
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
        delta = -JF(u,v)\F(u,v);        % remove ";" to see how delta decreases
        its = its + 1;
    end
   
    sol = sol+delta;
    u = sol(1:my);                      % it will became the next u_Eprevious
    v = sol(my+1:2*my);
    u_Eprevious = u;

end 

plot(u,y,v,y)
xlabel('Velocity');
ylabel('y');
text(2,2.5,'u', 'Fontsize',13)
text(0.9,3.5,'v','Fontsize',13)
legend('u', 'v');
legend('Location', 'SouthEast');
legend('Fontsize',12)
title('VELOCITY OF BOUNDING SURFACE (2D) WITH EXTRENAL FLUX (U,0)')
axis([0 3 0 10]);

toc

%% Check results
% Shape factor H = delta/theta with:

% delta = integral from 0 to +infinity of (1 - u(y)/U) dy
delta = trapz(y,1-u/U);

% theta = integral from 0 to +infinity of (1 - u(y)/U) * u(y)/U dy
theta = trapz(y,(1-u/U).*(u/U));

H = delta/theta;
fprintf('Shape factor is H = %f\n', H)
disp('Blasius boundary layer value of shape factor is H_boundary_layer = 2.5916')

H_Blasius = 2.5916;
if (abs(H - H_Blasius) < 0.05)
    disp('H - H_boundary_layer < 5%')
end
