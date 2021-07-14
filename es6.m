% ESERCIZIO 6

% BLASIUS EQUATION:
% NUMERICAL SOLUTION OF BOUNDARY LAYER OF A BOUNDING SURFACE (2D)
% During the calculation of Blasius equation, in momentum equation 
% we obtain:
% U/nu*g(x)*g'(x) = - f'"(eta)/(f(eta)*f"(eta)) and this is satisfied IFF both
% are equal to a positive constant called alpha.
% Let alpha be equal to 1/2, then:
% eta = y * sqrt(U/(nu*xb));
% f'" + 1/2 f*f" = 0 (from which we get the diffeq problem below through the 
% sostitution f' = u) with its boundary conditions

%%
% Solve the diffeq problem
% { fu' + 2u" = 0
% { f' - u = 0
% { f(0) = 0; f'(inf) = 1
% { u(0) = 0; u(inf) = 1
% using discretisation with finite difference

clear all
close all

tic
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

%__________________________ similarity variable eta
eta = y * sqrt(U/(nu*xb));


%_____________________________________ F1 and U1 

% denominator of subdiagonal and superdiagonal terms of finite
% differece matrix
d = 1./(h(1:m-2)+h(2:m-1));

% finite difference matrix for first derivative related to f
F1 = spdiags([[-d;0;0],[0;0;d]],[-1,1],m,m);

% finite difference matrix for first derivative related to u
U1 = spdiags([[-d;0;0],[0;0;d]],[-1,1],m,m); % remark: first and last lines are zeros lines.

%_____________________________________ F2 and U2 
% denominator of subdiagonal terms
d1 = 1./(h(1:m-2).*(h(1:m-2)+h(2:m-1)));

% denominator of superdiagonal terms 
d2 = 1./(h(2:m-1).*(h(1:m-2)+h(2:m-1)));

% finite difference matrix for second derivative related to u
% remark: it's common to put 1 position (1,1) and (m,m) to sodisfy boundary
% conditions
U2 = spdiags([[2*d1;0;0],[0;-2*d1-2*d2;0],[0;0;2*d2]],[-1,0,1],m,m);

%_____________________________________ known term
b = [zeros(2*m,1)];

%_____________________________________ boundary conditions (Dirichlet + Neumann)
% The discretisation of this system is:
% { f.*U1*u + 2*U2*u = 0
% { F1*f - u = 0
% and we have to sodysfy boundary conditions in the following rows:
% { f1.*U1*u1 + 2*U2*u1 = 0         (line 1)
% { fm.*U1*um + 2*U2*um = 0         (line m)
% { F1*f1 - u1 = 0                  (line m+1)
% { F1*fm - um = 0                  (line 2m)

% Leaving U1 the way it is, lines 1 and m become:
% { 2*U2*u1 = 0     ---> { u1 = 0         
% { 2*U2*um = 0     ---> { um = 0 so we need to change b(m) from 0 to 1   
% so we need to assign U2(1,1) = 1/2; U2(m,m) = 1/2. (It would be right also 
% U2(1,1) = 1 ;U2(m,m) = 1)
U2(1,1) = 1/2;
U2(m,m) = 1/2;
b(m) = 1;

% For Neumann condition we use the following approximation of first
% derivative 
% f'm = (fm - fm-1)/(h(end));  
% We know the value of f'm, so f'm = (fm - fm-1)/(h(end)) ---> 
% 1 = (fm - fm-1)/(h(end)) 
% Putting this identity in f' - u = 0, we obtain:
%  (fm - fm-1)/(h(end)) - um = 0 ---> (fm - fm-1)/(h(end)) = um = 1
%  and this is correct.
% In the last row of F1 we have to put 
F1(m,m-1) = -1./h(end);     
F1(m,m) = 1./h(end);      

% For last boundary condition we have
% { F1*f1 - u1 = 0 ---> { F1*f1 = 0 so we have to impose F1(1,1) = 1 to
% satisfy f1 = 0 
F1(1,1) = 1;


% Function
F = @(f,u) [f.*(U1*u) + 2*(U2*u) - b(1:m);...
    F1*f - u - b(m+1:2*m)];

% Jacobian 
JF = @(f,u) [diag(U1*u), f.*U1 + 2*U2;...
    F1, -eye(m)];

%_____________________________________ Newton
tol = 1e-8;
% initial guess should be near to the solution, so let u sodisfies
% boundary conditions (for f we don't have solution in m)
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

%_______________________________________ plot

plot(f,eta,u,eta)
text(1.2,6,'f^{\prime}(\eta)', 'Fontsize',13)
text(2.8,4.2,'f(\eta)','Fontsize',13)
xlabel('Velocity');
ylabel('Eta');
legend('f(eta)', 'u(eta)','Location', 'SouthEast');
legend('Fontsize',12)



%% Calculus of u(2,y) and v(2,y) where 2 is x_final = xb
%  Known facts:
%  u(x,y) = U * f'(eta);
%  v(x,y) = U * g'(x) * (eta * f'(eta) - f(eta));
%  where f'(eta) = u(eta) and
%  g(x) = sqrt(nu*x/U) ---> g'(x) = 1/2 *sqrt(nu/U) * 1/sqrt(x)

u_2 = U * u;
v_2 = U * 1/2 *sqrt(nu/U) * 1/sqrt(xb) * (eta.*u - f);

%_______________________________________ plot
figure
plot(u_2, y, v_2, y);
text(2,2.5,'u', 'Fontsize',13)
text(0.9,3.5,'v','Fontsize',13)
xlabel('Velocity');
ylabel('y');
legend('u', 'v');
legend('Location', 'SouthEast');
legend('Fontsize',12)
axis([0 3 0 10]);
toc

%% Check results
% Shape factor H = delta/theta with:

% delta = integral from 0 to +infinity of (1 - u(y)/U) dy
delta = trapz(y,1-u_2/U);

% theta = integral from 0 to +infinity of (1 - u(y)/U) * u(y)/U dy
theta = trapz(y,(1-u_2/U).*(u_2/U));

H = delta/theta;
fprintf('Shape factor is H = %f\n', H)
disp('Blasius boundary layer value of shape factor is H_Blasius = 2.5916')

H_Blasius = 2.5916;
if (abs(H - H_Blasius) < 0.05)
    disp('H - H_blasius < 5%')
end