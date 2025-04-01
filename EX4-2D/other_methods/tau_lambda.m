%% Riesz tau-preconditioner
function y=tau_lambda(a1,k)
C2=[a1(3:k);0;0]; % the first column of c's hankel matrix
C2=a1-C2;  % the first comun of tau matrix of Toeplitz A
s=(1:1:k)';
ss=sqrt(2/(k+1))*sin(pi/(k+1)*s);  % the first column of sine transform matrix
ss=1./ss; 
LAM1=sqrt(2/(k+1))*dst(C2); % S multiply the first of Tau matrix
y=ss.*LAM1;  % the eigenvalue of tau matrix 
end