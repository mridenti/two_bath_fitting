function [beta, Res, beta_exp, Res_exp, TT, X_fit, Y_fit] =  ...
    fit_OH(X, Y, Sy, X_all, Y_all, Sy_all)

% Two-exponential model was linearized applying a log function  

close all

global T kB;

global xi;

xi = 30;

% Physical constants definition
kB = 8.61733034E-05;

is_2T = 1; % set 0 if T1 should be kept constant / set 1 otherwise

% This number is model dependent; see comments below
m0 = 6; 
% 5-n_g for the most general exothermic case,
% 6-n_g for the general 2T model
% 9-n_g for general 2T model with 1-g
% 3-n_g for fixed width model or fixed power law, 
% 4-n_g for 2-T model or power law

% Set here the number of Gaussians to form ripples. Gaussian term order n_g 
n_g = 0;
% This formula sets the number of parameter of the model
n_g_param = m0 + (n_g-1)*3;

X = X(:);
Y = Y(:);
Sy = Sy(:);

X_all = X_all(:);
Y_all = Y_all(:);
Sy_all = Sy_all(:);

x = X;
y = Y;
sy = Sy;
n_data = length(x);
n_data_all = length(X_all);
d_data = n_data_all - n_data;

% Initial estimate of the temperature T1
x_l = X(1:6);
y_l = log(Y(1:6));
%sy_l = Sy(1:6);

[ Aj_lin ] = mmq([0 1], x_l, y_l);
[ Aj_lin ] = mmq([0 1], x_l, y_l, sqrt(sum(Aj_lin.Dif.^2)/Aj_lin.ngl));
a_l = Aj_lin.A(2);
sa_l = Aj_lin.SA(2);
b_l = Aj_lin.A(1);

% Initial temperature and error estimates
T(1) = 1/(kB*abs(a_l));
T(2) = kB*T(1)^2*sa_l;

% Initial estimate of the temperature T2
x_u = X(7:end);
y_u = log(Y(7:end));
%sy_l = Sy(1:6);

[ Aj_lin ] = mmq([0 1], x_u, y_u);
[ Aj_lin ] = mmq([0 1], x_u, y_u, sqrt(sum(Aj_lin.Dif.^2)/Aj_lin.ngl));
a_u = Aj_lin.A(2);

% Initial temperature and error estimates
T_u = 1/(kB*abs(a_u));

model = @(beta, x0) fTTe_xi(beta, x0);
model_T = @(beta, x0) fTTe_T(beta, x0);
model_2exp = @(beta, x) beta(1)*exp(-x/beta(3)) + beta(2)*exp(-x/beta(4));
log_model_2exp = @(beta, x) log(beta(1)*exp(-x/beta(3)) + beta(2)*exp(-x/beta(4)));
%beta0 = [1.5018e+73 6.1416e+70 5.4822e+00 6.1213e-01 4.7173e+04 4.9790e+00 2.4014e-01 9.6783e-01 kB*2.6245e+02]; % general 2-T model with 1-g
%beta0 = [exp(b_l) 2.0E-02*exp(b_l) 5.5e+00 3.0e-01 8.0e03 5.0 0.26 1.32
%kB*T(1)]; % general 2-T model with 1-g
beta0 = [exp(b_l) 4.2 0.06 kB*T_u kB*T(1)]; % general 2-T model
beta0_exp = [exp(b_l) 10 kB*T(1) kB*T_u];
%beta0 = [1.4E64 17.4 5.3 2.6e-02]; % power-law model
%beta0 = [0.003 5e65 1.4 2.6e-02]; % 2_T model
%beta0 = [ 0.8e+67   8.00e-03   4.0e+00   5.0e-01 2.5e-02 ]; % 1-g model
%beta0 = [ 4.8163e+66  2.1083e-03 4.1e+00 3.0929e-01 1.3758e-03 5.0e+00 8.7510e-01  2.5073e-02 ]; % 2-g model
%beta0 = [ exp(b_l)  5.0e-03   4.20000e+00   3.0e-01   2.0e-03  4.9e+00   1.0e-01  2.0e-03  5.6e+00   3.0e-01  kB*T(1)]; %3-g model
%beta0 = [ 1.0e+67   5.0e-03   3.0e-03   3.0e-03   4.10e+00   5.0e+00   6.0e+00 3.0e-01 2.5e-02]; % 3-gm model
x_fit = linspace(0.99*min(x), 1.01*max(x), 100);
y_fit = model(beta0, x_fit);
y_fit_exp = model_2exp(beta0_exp, x_fit);
%ngl_fit = length(x) - 5;
%ngl_fit_exp = length(x) - 4;
% plot the initial guess 
figure(1)
semilogy(x_fit,y_fit, '-r');
hold on;
semilogy(x_fit,y_fit_exp, '-g');
hold on;
% plot selected data points
errorbar(x,y,sy,'.');
hold on;
semilogy(x,exp(x*a_l+b_l),'-b');
hold on
% plot all data points, if some of them have been excluded
if d_data > 0
    errorbar(X_all(n_data+1:end),Y_all(n_data+1:end),Sy_all(n_data+1:end),'xr');
end
hold off;
maxiter = 30;
% Fit two bath
[beta, Res] = fitGM(x, y, 1, model, beta0, maxiter);
index_1 = (Res.x < 4.4); 
index_2 = (Res.x > 4.4); 
index_3 = (Res.x > 4.95);
ngl_fit_1 = length(Res.Res(index_1)) - 5;
ngl_fit_2 = length(Res.Res(index_2)) - 5;
ngl_fit_3 = length(Res.Res(index_3)) - 5;
corr_sy_1 = sqrt(sum(Res.Res(index_1).^2)/(ngl_fit_1));
corr_sy_2 = sqrt(sum(Res.Res(index_2).^2)/(ngl_fit_2));
corr_sy_3 = sqrt(sum(Res.Res(index_3).^2)/(ngl_fit_3));
sy_corr = sy;
sy_corr(index_1) = corr_sy_1;
sy_corr(index_2) = corr_sy_2;
sy_corr(index_3) = corr_sy_3;
[beta, Res] = fitGM(x, y, sy_corr, model, beta0, maxiter);
% Fit two exponential
[beta_exp, Res_exp] = fitGM(x, log(y), 1 , log_model_2exp, beta0_exp, maxiter);
%corr_sy_exp = sqrt(Res_exp.qui2/ngl_fit_exp);
%sy_corr_exp = corr_sy_exp*sy;
sy_corr_exp = sy_corr;
[beta_exp, Res_exp] = fitGM(x, log(y), sy_corr_exp./y, log_model_2exp, beta0_exp, maxiter);
T(1) = beta(n_g_param)/kB;
sy = sy_corr;
% Plot the fitted curve with data - semilogy scale
figure(2)
y_fit = model(beta,x_fit);
y_fit_exp = model_2exp(beta_exp, x_fit);
semilogy(x_fit,y_fit, '-r');
hold on;
semilogy(x,y,'o');
hold off;
TT = T;
% If the least squares does not converge try another fit with prescribed
% T1. This should be tried only if is_2T = 1
if (Res.iter == maxiter && is_2T ~= 1)
    fprintf('Fixando a temperatura em T = %f K \n',T(1));
    beta0_T = beta(1:(n_g_param-1));
    [beta_T, Res] = fitGM(x, y, sy, model_T, beta0_T, maxiter);
    beta(1:(n_g_param-1)) = beta_T;
    y_fit = model(beta,x_fit);
end
% plot fitted curve with data - linear scale
figure(3)
plot(x_fit,y_fit, '-r');
hold on;
plot(x_fit,y_fit_exp, '-g');
hold on;
errorbar(x,y,sy,'o');
hold off;
% plot fitted curve with data - semilogy scale
figure(4)
semilogy(x_fit,y_fit, '-r');
hold on;
semilogy(x_fit,y_fit_exp, '-g');
hold on;
errorbar(x,y,sy,'o');
hold off;
figure(5)
errorbar(x, Res.Res, ones(1,length(x)), '.')
hold on;
plot(x, zeros(1, length(x)));
hold off
% goes back to figure 1 and plot fitted curve 
figure(1)
hold on;
semilogy(x_fit,y_fit, '-b');
hold on;
semilogy(x_fit,y_fit_exp, '-g');
hold off;
X_fit = x_fit;
Y_fit = y_fit;
end

% Superposition 2-T model
function value = f0(beta,x)

value = beta(2)*exp(-x/beta(4)) + beta(1)*exp(-x/beta(3));

end

% Max entropy 2-T model with one gaussian ripple
function value = fTTe_c(beta, x)

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = (beta(1)+g0(x(i),beta, beta(9)))*...
        gT(x(i), beta(5), beta(6), beta(7), beta(9))*gTe(x(i), beta(5), beta(6),beta(7),beta(8));
end

end

% Max entropy 2-T model
function value = fTTe(beta, x)

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = beta(1)*gT(x(i), beta(2), beta(3), beta(4), beta(6))*gTe(x(i), beta(2), beta(3),beta(4),beta(5));
end

end

% xi is fixed
function value = fTTe_xi(beta, x)

global xi;

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = beta(1)*gT(x(i), xi, beta(2), beta(3), beta(5))*gTe(x(i), xi, beta(2),beta(3),beta(4));
end

end

% Power-law model
function value = fe(beta, x)

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = (beta(1) + ge(x(i), beta(2), beta(3), beta(4)))*exp(-x(i)/beta(4));
end

end

% One-gaussian model
function value = f1(beta, x)

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = (beta(1) + g1(x(i), beta(2), beta(3), beta(4), beta(5)))*exp(-x(i)/beta(5));
end

end
 
% Two gaussian model
function value = f2(beta, x)

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = (beta(1) + g2(x(i), beta(2), beta(3), beta(4), beta(5), beta(6), beta(7),beta(8)))*exp(-x(i)/beta(8));
end

end

% Three gaussian model - correction to prevent convergence to unphysical
% mean energies
function value = f3(beta, x)

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = (beta(1) + g3(x(i), beta(2), beta(3), beta(4), beta(5), beta(6), beta(7),beta(8), beta(9), beta(10), beta(11)))*exp(-x(i)/beta(11));
end

end

% Three gaussian model without correction
function value = f3m(beta, x)

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = (beta(1) + g3m(x(i), beta(2), beta(3), beta(4), beta(5), beta(6), beta(7),beta(8), beta(9)))*exp(-x(i)/beta(9));
end

end

% Max entropy 2-T model with one gaussian ripple - T1 is prescribed
function value = fTTe_cT(beta, x)

global T kB;

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = (beta(1)+g0(x(i),beta, kB*T(1)))*...
        gT(x(i), beta(5), beta(6), beta(7), kB*T(1))*gTe(x(i), beta(5), beta(6),beta(7),beta(8));
end

end

% Max entropy 2-T model - T1 is prescribed
function value = fTTe_T(beta, x)

value = zeros(length(x),1);

global T kB;

for i = 1:length(x)
    value(i) = beta(1)*gT(x(i), beta(2), beta(3), beta(4), kB*T(1))*gTe(x(i), beta(2), beta(3),beta(4),beta(5));
end

end

% Superposition 2-T model - T1 is prescribed 
function value = f0_T(beta,x)

global T kB;

value = beta(2)*exp(-x/(kB*T(1))) + beta(1)*exp(-x/beta(3));

end

% Power-law model - T1 is prescribed
function value = fe_T(beta, x)

global T kB;

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = (beta(1) + ge(x(i), beta(2), beta(3), kB*T(1)))*exp(-x(i)/(kB*T(1)));
end

end

% One-gaussian model - T1 is prescribed
function value = f1_T(beta, x)

global T kB;

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = (beta(1) + g1(x(i), beta(2), beta(3), beta(4), kB*T(1)))*exp(-x(i)/(kB*T(1)));
end

end

% Two-gaussian model - T1 is prescribed
function value = f2_T(beta, x)

global T kB;

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = (beta(1) + g2(x(i), beta(2), beta(3), beta(4), beta(5), beta(6), beta(7), kB*T(1)))*exp(-x(i)/(kB*T(1)));
end

end

% Three-gaussian model - T1 is prescribed. Includes algorithm to prevent
% unphysical mean energies 
function value = f3_T(beta, x)

global T kB;

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = (beta(1) + g3(x(i), beta(2), beta(3), beta(4), beta(5), beta(6), beta(7), beta(8), beta(9), beta(10), kB*T(1)))*exp(-x(i)/(kB*T(1)));
end

end

% Three-gaussian model - T1 is prescribed
function value = f3m_T(beta, x)

global T kB;

value = zeros(length(x),1);

for i = 1:length(x)
    value(i) = (beta(1) + g3m(x(i), beta(2), beta(3), beta(4), beta(5), beta(6), beta(7), beta(8), kB*T(1)))*exp(-x(i)/(kB*T(1)));
end

end

function value = gT(x, b2, b3, b4, b6)

y = @(xx) 1./((1+b2*normcdf(xx,b3,b4))*b6);
value = exp(-integral(y,0,x));

end

function value = gTe(x, b2, b3, b4, b5)

y = @(xx) 1./((1+(1./(b2*normcdf(xx,b3,b4))))*b5);
value = exp(-integral(y,0,x));

end

function value = ge(x, b2, b3, b4)

y = @(xx) b2*(xx.^(-b3)).*exp(xx/b4);
value = integral(y,0,x);


end

function value = g0(x, beta, T)

y = @(xx) beta(2)*exp(-(xx-beta(3)).^2/(2*beta(4)^2)).*exp(xx/T);
value = integral(y,0,x);

end

function value = g1(x, b2, b3, b4, b5)

y = @(xx) (b2*exp(-(xx-b3).^2/(2*b4^2))).*exp(xx/b5);
value = integral(y,0,x);

end

function value = g2(x, b2, b3, b4, b5, b6, b7, b8)

y = @(xx) (b2*exp(-(xx-b3).^2/(2*b4^2))+b5*exp(-(xx-b6).^2/(2*b7^2))).*exp(xx/b8);
value = integral(y,0,x);

end

% 3-gaussian integral. Algorithm to prevent convergence to unphysical
% mean energies  
function value = g3(x, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11)

y = @(xx) (b2*exp(-(xx-b3).^2/(2*b4^2))+b5*exp(-(xx-b6).^2/(2*b7^2))+b8*exp(-(xx-b9).^2/(2*b10^2))).*exp(xx/b11);
value = integral(y,0.0,x);
e_l = 4.01;
if(b3  < e_l) 
    value = (1+exp(abs(b3-e_l)))*value;
end
if(b6  < e_l) 
    value = (1+exp(abs(b6-e_l)))*value;
end
if(b9  < e_l)
    value = (1+exp(abs(b9-e_l)))*value;
end

end

function value = g3m(x, b2, b3, b4, b5, b6, b7, b8, b9)

y = @(xx) (b2*exp(-(xx-b5).^2/(2*b8^2))+b3*exp(-(xx-b6).^2/(2*b8^2))+b4*exp(-(xx-b7).^2/(2*b8^2))).*exp(xx/b9);
value = integral(y,0,x);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Levenberg-Gauss-Marquardt algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ beta, Res ] = fitGM( X, y, s, model, beta0, maxiter)
% [ beta, Res ] = mmqGM( X, y, s, model, beta0 )
% Argumentos de entrada:  
%   beta0: parâmetros iniciais
%
%
%
%  inspirado na rotina do nlinfit do toolbox do matlab
%  Primeira versao: Otaviano (14-9-01)
%  Ultima modificacao: Marco (05-05-09)

n = length(y);
%Transforma os vetores de dados em vetores coluna
y = y(:);
X = X(:);

p = length(beta0);
%Transforma os vetores de parâmetros iniciais em vetores coluna
beta0 = beta0(:);
if length(s)==1
    S = s*ones(n,1);
else
    S = s;
end

J = zeros(n,p);
beta = beta0;
betanew = beta + 1;
%maxiter = 40; % ==============
iter = 0;
betatol = 1.0E-2;
rtol = 1.0E-4;
qui2 = 1;
qui2old = qui2;
S = S(:);
v = eye(p);
convergiu = 0;
fatorGM = 1;
% disp( 'new fit' )
usandoGM = 0;
while and( convergiu<2, iter < maxiter )
    
    if iter > 0
        beta = betanew;
        yfit = yfitnew;
        r = rnew;
        qui2old = qui2;
    else
        yfit = feval(model,beta,X);
        r = ( y(:) - yfit(:) ) ./ S;
        qui2old = r'*r;
    end
   iter = iter + 1;
   J = zeros(size(y,1), p);
   for k = 1:p,
       delta = zeros( size( beta) );
       delta(k) = sqrt(eps)*beta(k);
        if abs(delta(k)) < 100*eps
           delta(k) = 100*eps;
        end
        if (abs(delta(k)) < 0.001*abs(beta(k)))
            delta(k) = 0.001*beta(k);
        end
       betaplus = beta + delta;
       betaminus = beta - delta;
       yplus = feval(model, betaplus , X);
       yminus = feval(model, betaminus, X);
       J(:,k) = ( ( yplus(:) - yminus(:) ) ./ S) / (delta(k));
       %J(:,k) = ( ( yplus(:) - yminus(:) ) ./ s ) / (betaplus(k)-betaminus(k));
   end
   JJ = J'*J;

   JJplus = JJ + (diag(diag(JJ)))*fatorGM;
   step = inverta(JJplus)*(J'*r);
   

   betanew = beta + step;
   yfitnew = feval(model,betanew,X);
   rnew = (y(:) - yfitnew(:))./S;
   
   qui2 = rnew'*rnew;
   if ( qui2 > qui2old*(1+1E-4))
       convergiu = 0;
       while and( qui2>qui2old*(1+1E-4), fatorGM<=1E2 )
           fatorGM = fatorGM * sqrt(10);
           JJplus = JJ + diag(diag(JJ))*fatorGM;
           step = inverta(JJplus)*(J'*r);
           betanew = beta + step;
           yfitnew = feval(model,betanew,X);
           rnew = (y(:) - yfitnew(:))./S;
           qui2 = rnew'*rnew;
       end
       fatorGM = 1;
       %       disp( sprintf( '%3d', iter) );
   else
       fatorGM = fatorGM / 10;
       if abs( qui2old - qui2 ) < 1E-5*(qui2+qui2old) %originalmente 1E-5*
           v = inverta( JJ )*eye(p);
           if max( abs( step )./sqrt( diag(v)*(qui2/(n-p)) ) )<1E-1 
               convergiu = convergiu +1;
               fatorGM = 1E-6;
           end
       end
   end
end
if convergiu==0, v = inverta( JJ )*eye(p); end
if iter == maxiter
   disp('mmqGM nao convergiu. Retornando os resultados da ultima iteracao.');
end
Res.A = beta;
Res.SA = sqrt( diag(v) );
Res.VA = v;
Res.qui2 = qui2;
Res.pqui2 = 100*pqui2(qui2, n-p);
Res.ngl = n-p;
Res.F = yfit;
Res.iter = iter;
Res.si = S;
Res.x = X;
Res.y = y;
Res.Dif = y(:) - yfit(:);
Res.Res = r;
end

function [ invM ] = inverta( M )
% [ invM ] = inverta( M )
%
% usar no lugar da função inv

% (c) Zwinglio Guimaraes-Filho (2001-2002)

% para verificar se há elementos na diagonal principal menores do que a
% precisão numérica (eps):
%
%  if min( sqrt( diag(M) ) ) < eps
%  M
%  diag(M)
%  error('ops')
%  end

dM = 1./( sqrt( diag(M) ) );
MdM = dM * dM';

invMuni = inv( MdM .* M );
invM = MdM .* invMuni;

end

function resultado = pqui2(chi2,N)

if (chi2<=0 || ~isscalar(chi2))
    disp('Chi2 deve ser estritamente positivo e escalar.');
    return;
end

if (N<1 || ~isscalar(N))
    disp('Número de graus de liberdade N deve ser inteiro e escalar.');
    return;
end

    try
        resultado = 1 - chi2cdf(chi2 , N );
    catch
        F = @(x) ((0.5*x^2)^(0.5*N-1))*exp(-0.5*x^2)/gamma(0.5*N);
        P = 1 - quadl(F, 0, 10*N);
    end
end

function p = chi2cdf(x,v)
%CHI2CDF Chi-square cumulative distribution function.
%   P = CHI2CDF(X,V) returns the chi-square cumulative distribution
%   function with V degrees of freedom at the values in X.
%   The chi-square density function with V degrees of freedom,
%   is the same as a gamma density function with parameters V/2 and 2.
%
%   The size of P is the common size of X and V. A scalar input   
%   functions as a constant matrix of the same size as the other input.    

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.4.

if   nargin < 2, 
    error('Requires two input arguments.');
end

[errorcode, x, v] = distchck(2,x,v);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end
    
% Call the gamma distribution function. 
p = gamcdf(x,v/2,2);

% Return NaN if the degrees of freedom is not a positive integer.
k = find(v < 0  |  round(v) ~= v);
if any(k)
   tmp = NaN;
    p(k) = tmp(ones(size(k)));
end

end