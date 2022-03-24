% Eigenvalue Search for the Equations of Fiber Spinning
%
%		      Md Fayaz Ahamed
%             February 2021
% ------------------------------------------------------------
clc
tic
%diary on
% 16 digit floating point
format long e
D=21;
%parameters
beta=1.0;ita=1.0;F=1.0;
% N	- number of collocation points -1 (0,..,N)
N = input('Please enter the number of collacation point N you want:');
% collocation points
X = 0:N;
Xa = cos(pi*X/N);
Xb = 0.5*(1.0+Xa(1:N));
% ------------------------------------------------------------

% Chebyshev polynomials
T = zeros(N+1);
DT = zeros(N+1);
%
T(:,1) = 1;
T(:,2) = Xa';
DT(:,1) = 0;
DT(:,2) = 1;
%
for i = 3:N+1
	T(:,i) = 2*Xa'.*T(:,i-1)-T(:,i-2);
	DT(:,i) = 2*T(:,i-1)+2*Xa'.*DT(:,i-1)-DT(:,i-2);
end
T1=T(N+1,:);

% ------------------------------------------------------------
%
% Aim: Find the roots of det(A-B*lambda) = 0
%
Bmatrix = zeros(3*(N+1));
%
% Elements for b-equation
Bmatrix(1:N,1:(N+1)) = T(1:N,:);
Bmatrix((N+1):(N+N),(N+2):(N+N+2)) = T(1:N,:);
%
V=exp((F/4.0*ita)*Xb(1:N));
VD=(F/4.0*ita)*(exp((F/4.0*ita)*Xb(1:N)));
A=exp(-(F/4.0*ita)*Xb(1:N));
AD=(-F/4.0*ita)*(exp(-(F/4.0*ita)*Xb(1:N)));
B=exp(-(F/4.0*ita)*Xb(1:N)) + ((4.0*beta)/F)*(exp(-(F/4.0*ita)*Xb(1:N))).*((exp(-(F/4.0*ita)*Xb(1:N))) - 1.0);
BD=-(F/4.0*ita)*exp(-(F/4.0*ita)*Xb(1:N))-(2.0*beta/ita)*exp(-(F/2.0*ita)*Xb(1:N))+(beta/ita)*exp(-(F/4.0*ita)*Xb(1:N));

%
Vbar=diag(V);
Vbarprime=diag(VD);
Abar=diag(A);
Abarprime=diag(AD);
Bbar=diag(B);
Bbarprime=diag(BD);

% Entiresfor A-matrix
P= 2.0*Vbar*DT(1:N,:)+Vbarprime*T(1:N,:);
Q=Abarprime*T(1:N,:)+2.0*Abar*DT(1:N,:);
R=Bbarprime*T(1:N,:)+2.0*Bbar*DT(1:N,:);
S=2.0*Bbar*Vbarprime-2.0*Bbar*Vbar+(beta/ita)*Abar;
L =3.0*Abar*Vbarprime+1.0-2.0*Vbar*Abarprime-(F/ita);
K =3.0*Abar*Bbar + Abar*Bbarprime + 2.0*Bbar*Abarprime;
ST=S*T(1:N,:);
LT=L*T(1:N,:);
KT=K*T(1:N,:);
%
Amatrix = zeros(3*(N+1));
%
% Equations for determining elements of Amatrix
Amatrix(1:N,1:(N+1)) = -P;
Amatrix(1:N,(N+N+3):(N+N+N+3))= -Q;
Amatrix((N+1):(N+N),1:(N+1))=-(beta/ita)*T(1:N,:);
Amatrix((N+1):(N+N),(N+2):(N+N+2)) = -P; 
Amatrix((N+1):(N+N),(N+N+3):(N+N+N+3))= -R;
Amatrix((2*N+1):(3*N),1:(N+1))=ST;
Amatrix((2*N+1):(3*N),(N+2):(2*N+2)) = LT; 
Amatrix((2*N+1):(3*N),(2*N+3):(3*N+3))=KT;
Amatrix((3*N+1),1:(N+1))=T1; 
Amatrix((3*N+2),(N+2):(2*N+2))=T1; 
Amatrix((3*N+3),(2*N+3):(3*N+3))=T1; 
Amatrix;
% ------------------------------------------------------------
%
% Eigenvalue calculation - filtering
%
eigenvalues = eig(Amatrix,Bmatrix)
plot(eigenvalues,'o')
xlabel('Real')
ylabel('Imaginary')
%
toc
% index = find(finite(eigenvalues));
% result = eigenvalues(index);
% [realparts,ind] = sort(real(result));
% while any(realparts>  0.5*nu*exp(3*nu)/(exp(nu)-1))
% 	ind = ind(1:size(ind)-1);
% 	realparts = realparts(1:size(ind));
% end
% %
% % Array of eigenvalues
% complex_eigenvalues = result(ind)
%
%diary off
%
% ------------------------------------------------------------
% End of code


