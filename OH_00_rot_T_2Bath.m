function [beta, Res, Yerror, error95CI] = OH_00_rot_T_2Bath

close all;

% plot experimental data

figure(1);
file_id_entrada = fopen('entrada.dat','r');
fgetl(file_id_entrada);
ExpData = fscanf(file_id_entrada,'%g',[1 inf]);
ExpData = .1*ExpData';
XData = ExpData;
file_id_saida = fopen('saida.dat','r');
fgetl(file_id_saida);
ExpData = fscanf(file_id_saida,'%g',[1 inf]);
ExpData = .1*ExpData';
YData = ExpData;
% select OH SPS - (1,3)/(0,2) 
xlim_sup = 312.0; % 312.0
xlim_inf = 306.0; % 306.0
YData = YData(XData<xlim_sup & XData>xlim_inf);
XData = XData(XData<xlim_sup & XData>xlim_inf);

% XData correction
XData = XData - (XData(YData == max(YData(XData < 307.96 & XData >307.8)))) + 307.93375;
YData = YData - mean(YData(XData<306.3));
%YData = YData - mean(YData(XData<382));
%YData = YData - min(YData);
YData = YData/max(YData);
plot(XData,YData);
axis([306 312 0 1]);
hold on;

global N_J;
global h_p c kb e_0 e m_e;
global cte;
global delta_lambda_inst;
global OH_En_X1_0f OH_En_X2_0f;
global OH_En_A1_0 OH_En_A2_0 OH_En_X1_0e  OH_En_X2_0e ... 
    OH_DEn_R1_00 OH_DEn_R2_00 OH_DEn_P1_00 OH_DEn_P2_00 ...
    OH_DEn_Q1_00 OH_DEn_Q2_00 ;
global A_OH_00 V_C_0 OH_a_X_0; 

delta_lambda_inst = 0.5e-02;
N_J = 100;

format long;

% physical constants
e_0 = 8.854187817E-12; %vacuum dielectric permitivity
m_e = 9.10938356E-31; %electron mass
c = 299792458; %light speed
e = 1.6021766208E-19;%electric charge
h_p = 6.626070040E-34;% Planck constant
kb = 1.38064852E-23;

% pre-factor in absorption coefficients
cte = e^2/(4*m_e*c*e_0);

% Spectroscopic constants and parameters
V_C_0 = 32440.60;

% Vibrational Einstein-coefficients 
A_OH_00 = 1.451E6;

% molecular constants

% OH A-X transition, v'-v'' (0 - 0) 
% Destombes (1977) - see function OH_Energy_00(J, p) in the last lines

% Lower State X 2Pi

% (nu = 0)

% upper state A 2Sigma

% (nu = 0)

% (0-0)
OH_En_A1_0 = zeros(N_J,1);
OH_En_A2_0 = zeros(N_J,1);
OH_En_X1_0e = zeros(N_J,1);
OH_En_X2_0e = zeros(N_J,1);
OH_En_X2_0f = zeros(N_J,1);
OH_En_X1_0f= zeros(N_J,1);
OH_DEn_R1_00 = zeros(N_J,1);
OH_DEn_R2_00 = zeros(N_J,1);
OH_DEn_P1_00 = zeros(N_J,1);
OH_DEn_P2_00 = zeros(N_J,1);
OH_DEn_Q1_00 = zeros(N_J,1);
OH_DEn_Q2_00 = zeros(N_J,1);

    for J=2:N_J+1
       
       % OH energies (solution of the secular equation)

        [OH_En_A2_0(J), OH_En_X2_0e(J), OH_En_X1_0e(J)] = OH_Energy_00(J-1-0.5,1);
        [OH_En_A1_0(J), OH_En_X2_0f(J), OH_En_X1_0f(J)] = OH_Energy_00(J-1-0.5,-1);
    end

    for J=2:N_J+1     
            
        % OH transition energies (0-0)
            
        OH_DEn_R1_00(J) = OH_En_A1_0(J) - OH_En_X1_0f(J-1);
        OH_DEn_R2_00(J) = OH_En_A2_0(J) - OH_En_X2_0e(J-1);                        
        OH_DEn_P1_00(J) = OH_En_A1_0(J-1) - OH_En_X1_0f(J);            
        OH_DEn_P2_00(J) = OH_En_A2_0(J-1) - OH_En_X2_0e(J); 
        OH_DEn_Q1_00(J) = OH_En_A1_0(J) - OH_En_X1_0e(J);
        OH_DEn_Q2_00(J) = OH_En_A2_0(J) - OH_En_X2_0f(J);
            
    end 

modelo = @(beta,lambda) I_syntetic(beta,lambda);
%beta0 = [1100 5000 2000 delta_lambda_inst 2000]; 
beta0 = [1118 4610 2100 2819 delta_lambda_inst 1.5E-22]; 
plot(XData,feval(modelo,beta0,XData),'r');
hold off;
[~, Res] = fitGM(XData, YData, 1, modelo, beta0, 30); 
Yerror = sqrt(Res.qui2/Res.ngl);
[beta, Res] = fitGM(XData, YData, Yerror, modelo, beta0, 30); 
[T1plus, T2plus, E0plus, T1minus, T2minus, E0minus] = ...
    evaluateUncertainty(XData, YData, Yerror, modelo, beta);
I = I_syntetic(beta,XData);
figure(2)
plot(XData,I);
axis([306 312 0 1]);
figure(3)
plot(XData,I,'-r');
hold on;
plot(XData,YData);
axis([306 312 0 1]);
hold off;
save output
error95CI = [T1plus, T2plus, E0plus, T1minus, T2minus, E0minus]; 
end

function [T1plus, T2plus, E0plus, T1minus, T2minus, E0minus] = ...
    evaluateUncertainty(XData, YData, Yerror, modelo, beta)
    beta_init = beta;
    ngl = length(XData);
    pqui2 = chi2cdf(ngl, ngl);
    ci = 0.68;
    while pqui2 < ci
        beta(1) = beta(1) + 0.001*beta(1);
        chi2_new = sum((feval(modelo, beta, XData)-YData).^2);
        pqui2 = chi2cdf(chi2_new/Yerror^2, ngl);
    end
    T1plus = beta(1);
    beta = beta_init;
    pqui2 = chi2cdf(ngl, ngl);
    while pqui2 < ci
        beta(1) = beta(1) - 0.001*beta(1);
        chi2_new = sum((feval(modelo, beta, XData)-YData).^2);
        pqui2 = chi2cdf(chi2_new/Yerror^2, ngl);
    end
    T1minus = beta(1);
    beta = beta_init;
    pqui2 = chi2cdf(ngl, ngl);
    while pqui2 < ci
        beta(2) = beta(2) + 0.001*beta(2);
        chi2_new = sum((feval(modelo, beta, XData)-YData).^2);
        pqui2 = chi2cdf(chi2_new/Yerror^2, ngl);
    end
    T2plus = beta(2);
    beta = beta_init;
    pqui2 = chi2cdf(ngl, ngl);
    while pqui2 < ci
        beta(2) = beta(2) - 0.001*beta(2);
        chi2_new = sum((feval(modelo, beta, XData)-YData).^2);
        pqui2 = chi2cdf(chi2_new/Yerror^2, ngl);
    end
    T2minus = beta(2);
    beta = beta_init;
    pqui2 = chi2cdf(ngl, ngl);
    while pqui2 < ci
        beta(3) = beta(3) + 0.001*beta(3);
        chi2_new = sum((feval(modelo, beta, XData)-YData).^2);
        pqui2 = chi2cdf(chi2_new/Yerror^2, ngl);
    end
    E0plus = beta(3);
    beta = beta_init;
    pqui2 = chi2cdf(ngl, ngl);
    while pqui2 < ci
        beta(3) = beta(3) - 0.001*beta(3);
        chi2_new = sum((feval(modelo, beta, XData)-YData).^2);
        pqui2 = chi2cdf(chi2_new/Yerror^2, ngl);
    end
    E0minus = beta(3);
end

function I = I_syntetic(beta, lambda)
global N_J;
global h_p c kb;
global delta_lambda_inst;
global OH_En_X1_0f OH_En_X2_0f;
global OH_En_A1_0 OH_En_A2_0 OH_En_X1_0e OH_En_X2_0e ... 
    OH_DEn_R1_00 OH_DEn_R2_00 OH_DEn_P1_00 OH_DEn_P2_00 ...
    OH_DEn_Q1_00 OH_DEn_Q2_00 ;
global A_OH_00 V_C_0 OH_a_X_0; 

lambda = lambda(:);
T1rotOH = abs(beta(1));
T2rotOH = abs(beta(2));
E0 = abs(beta(3));
%E0 = 1300;
s0 = abs(beta(4));
%s0 = 2850;
delta_lambda_inst = beta(5);
OH_ratio = 1;
sum = zeros(length(lambda),1);
    for J=2:N_J  
        Jm = J-0.5;
        gP_OH = 4*(2*(Jm-1)+1);
        gR_OH = 4*(2*(Jm+1)+1);
        gQ_OH = 4*(2*Jm+1);
        
        % OH frequencies and wavelengths
        
        % (0-0)        
        OH_nu_R1_00 = (OH_DEn_R1_00(J+1))*1E2*c;
        OH_nu_R2_00 = (OH_DEn_R2_00(J+1))*1E2*c;
        OH_nu_P1_00 = (OH_DEn_P1_00(J))*1E2*c;
        OH_nu_P2_00 = (OH_DEn_P2_00(J))*1E2*c;
        OH_nu_Q1_00 = (OH_DEn_Q1_00(J))*1E2*c;
        OH_nu_Q2_00 = (OH_DEn_Q2_00(J))*1E2*c;
     
        
        OH_lambda_P1_00 = 1E7/OH_DEn_P1_00(J);
        OH_lambda_R1_00 = 1E7/OH_DEn_R1_00(J+1);   
        OH_lambda_P2_00 = 1E7/OH_DEn_P2_00(J+1);
        OH_lambda_R2_00 = 1E7/OH_DEn_R2_00(J); 
                
        OH_lambda_Q1_00 = 1E7/OH_DEn_Q1_00(J);  
        OH_lambda_Q2_00 = 1E7/OH_DEn_Q2_00(J);
        
        
%        SJR1T_00 = SRDoublet1(Jm+1,1,N2_Y_B,N2_Y_C)/(4*(2*(Jm+1)+1));
%        SJR2T_00 = SRDoublet2(Jm+1,1,N2_Y_B,N2_Y_C)/(4*(2*(Jm+1)+1));
%        SJQ1T_00 = SQDoublet1(Jm,1,N2_Y_B,N2_Y_C)/(4*(2*(Jm)+1));
%        SJQ1T_00 = SQDoublet2(Jm,1,N2_Y_B,N2_Y_C)/(4*(2*(Jm)+1));
%        SJP1T_00 = SPDoublet1(Jm-1,1,N2_Y_B,N2_Y_C)/(4*(2*(Jm-1)+1));
%        SJP1T_00 = SPDoublet2(Jm-1,1,N2_Y_B,N2_Y_C)/(4*(2*(Jm-1)+1));

%         SJR1T_00 = 0.5*(Jm+2)/(4*(2*(Jm+1)+1));
%         SJR2T_00 = 0.5*(Jm+2)/(4*(2*(Jm+1)+1));
%         SJQ1T_00 = 0.5*(2*Jm+1)/(4*(2*Jm+1));
%         SJQ2T_00 = 0.5*(2*Jm+1)/(4*(2*Jm+1));
%         SJP1T_00 = 0.5*((Jm)-1)/(4*(2*(Jm-1)+1));
%         SJP2T_00 = 0.5*((Jm)-1)/(4*(2*(Jm-1)+1));

         SJR1T_00 = 0.5*(Jm+1.5)*(Jm+0.5)/(Jm+1)/(4*(2*(Jm+1)+1));
         SJR2T_00 = 0.5*(Jm-0.5)*(Jm+0.5)/(Jm+1)/(4*(2*(Jm+1)+1));
         SJQ1T_00 = 0.5*(Jm+0.5)^2*(2*Jm+1)/(Jm*(Jm+1))/(4*(2*Jm+1));
         SJQ2T_00 = 0.5*(Jm-0.5)*(Jm+1.5)*(2*Jm+1)/(Jm*(Jm+1))/(4*(2*Jm+1));
         SJP1T_00 = 0.5*(Jm-0.5)*(Jm+0.5)/Jm/(4*(2*(Jm-1)+1));
         SJP2T_00 = 0.5*(Jm+1.5)*(Jm+0.5)/Jm/(4*(2*(Jm-1)+1));
        
            % (nu 0-0) transition in OH R1, R2, R3, Q1, Q2 and Q3 branches
                  
        sum = sum + OH_ratio*A_OH_00*gQ_OH*SJQ2T_00*OH_nu_Q2_00*...
            exp(-c*h_p*(OH_En_A2_0(J)-V_C_0)*1E2*normcdf((E0 - OH_En_A2_0(J)+V_C_0)/s0)/kb/T1rotOH).*...
            exp(-c*h_p*(OH_En_A2_0(J)-V_C_0)*1E2*normcdf((OH_En_A2_0(J)-V_C_0 - E0)/s0)/kb/T2rotOH).*...
            exp(-(lambda-OH_lambda_Q2_00).^2*log(2)/delta_lambda_inst^2); 
        sum = sum + OH_ratio*A_OH_00*gR_OH*SJR2T_00*OH_nu_R2_00*...
            exp(-c*h_p*(OH_En_A2_0(J+1)-V_C_0)*1E2*normcdf((E0 - OH_En_A2_0(J+1)+V_C_0)/s0)/kb/T1rotOH).*...
            exp(-c*h_p*(OH_En_A2_0(J+1)-V_C_0)*1E2*normcdf((OH_En_A2_0(J+1)-V_C_0 - E0)/s0)/kb/T2rotOH).*...
            exp(-(lambda-OH_lambda_R2_00).^2*log(2)/delta_lambda_inst^2);  
        if (J>2)
            sum = sum + OH_ratio*A_OH_00*gR_OH*SJR1T_00*OH_nu_R1_00*...
                exp(-c*h_p*(OH_En_A1_0(J+1)-V_C_0)*1E2*normcdf((E0 - OH_En_A1_0(J+1)+V_C_0)/s0)/kb/T1rotOH).*...
                exp(-c*h_p*(OH_En_A1_0(J+1)-V_C_0)*1E2*normcdf((OH_En_A1_0(J+1)-V_C_0 - E0)/s0)/kb/T2rotOH).*...
                exp(-(lambda-OH_lambda_R1_00).^2*log(2)/delta_lambda_inst^2);
            sum = sum + OH_ratio*A_OH_00*gQ_OH*SJQ1T_00*OH_nu_Q1_00*...
                exp(-c*h_p*(OH_En_A1_0(J)-V_C_0)*1E2*normcdf((E0 - OH_En_A1_0(J)+V_C_0)/s0)/kb/T1rotOH).*...
                exp(-c*h_p*(OH_En_A1_0(J)-V_C_0)*1E2*normcdf((OH_En_A1_0(J)-V_C_0 - E0)/s0)/kb/T2rotOH).*...
                exp(-(lambda-OH_lambda_Q1_00).^2*log(2)/delta_lambda_inst^2);
            sum = sum + OH_ratio*A_OH_00*gP_OH*SJP1T_00*OH_nu_P1_00*...
                exp(-c*h_p*(OH_En_A1_0(J-1)-V_C_0)*1E2*normcdf((E0 - OH_En_A1_0(J-1)+V_C_0)/s0)/kb/T1rotOH).*...
                exp(-c*h_p*(OH_En_A1_0(J-1)-V_C_0)*1E2*normcdf((OH_En_A1_0(J-1)-V_C_0 - E0)/s0)/kb/T2rotOH).*...
                exp(-(lambda-OH_lambda_P1_00).^2*log(2)/delta_lambda_inst^2);
            sum = sum + OH_ratio*A_OH_00*gP_OH*SJP2T_00*OH_nu_P2_00*...
                exp(-c*h_p*(OH_En_A2_0(J-1)-V_C_0)*1E2*normcdf((E0 - OH_En_A2_0(J-1)+V_C_0)/s0)/kb/T1rotOH).*...
                exp(-c*h_p*(OH_En_A2_0(J-1)-V_C_0)*1E2*normcdf((OH_En_A2_0(J-1)-V_C_0 - E0)/s0)/kb/T2rotOH).*...
                exp(-(lambda-OH_lambda_P2_00).^2*log(2)/delta_lambda_inst^2);           
        end
        
    end
    I = beta(6)*sum;
end

function f = ff(x)
    if abs(x+1)<=1
        f = (x+2);
    elseif x < -2
        f = 0;
    else
        f= 1;
    end
   
end

function f = ff0(x)
    if x < 0
        f = 0;
    else
        f= 1;
    end
   
end

function [ beta, Res ] = fitGM( X, y, s, model, beta0, maxiter)

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
if min( sqrt( diag(M) ) ) < eps
M
diag(M)
error('ops')
end

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

% General Honl-London Formulas in the Born-Op approx. 

function u = uP(J,Y,L)
    u = sqrt(L^2*Y*(Y-4)+4*(J+0.5)^2)+L*(Y-2);
end

function u = uM(J,Y,L)
    u = sqrt(L^2*Y*(Y-4)+4*(J+0.5)^2)-L*(Y-2);
end

function C = CP(J,Y,L)
    C = 0.5*(uP(J,Y,L)^2+4*((J+0.5)^2-L^2));
end

function C = CM(J,Y,L)
    C = 0.5*(uM(J,Y,L)^2+4*((J+0.5)^2-L^2));
end

function S = SPDoublet1(J, L, Y1, Y2)

S = (J-L-0.5)*(J+L+0.5)*(uM(J-1,Y1,L)*uM(J-1,Y2,L)+...
    4*(J-L+0.5)*(J+L-0.5))^2/...
    (4*J*CM(J-1,Y1,L)*CM(J,Y2,L));

end

function S = SPDoublet2(J, L, Y1, Y2)

S = (J-L-0.5)*(J+L+0.5)*(uP(J-1,Y1,L)*uP(J,Y2,L)+...
    4*(J-L+0.5)*(J+L-0.5))^2/(4*J*CM(J-1,Y1,L)*CM(J,Y2,L));

end

function S = SRDoublet1(J, L, Y1, Y2)

S = (J-L+0.5)*(J+L+1.5)*(uM(J+1,Y1,L)*uM(J,Y2,L)+...
    4*(J-L+1.5)*(J+L+0.5))^2/(4*(J+1)*CM(J+1,Y1,L)*CM(J,Y2,L));

end

function S = SRDoublet2(J, L, Y1, Y2)

S = (J-L+0.5)*(J+L+1.5)*(uP(J+1,Y1,L)*uP(J,Y2,L)+...
    4*(J-L+1.5)*(J+L+0.5))^2/(4*(J+1)*CP(J+1,Y1,L)*CP(J,Y2,L));

end

function S = SQDoublet1(J, L, Y1, Y2)

S = (J+0.5)*((L+0.5)*uM(J,Y1,L)*uM(J,Y2,L)...
    +4*(L-0.5)*(J-L+0.5)*(J+L+0.5))^2/...
    (2*J*(J+1)*CM(J,Y1,L)*CM(J,Y2,L));

end

function S = SQDoublet2(J, L, Y1, Y2)

S = (J+0.5)*((L+0.5)*uP(J,Y1,L)*uP(J,Y2,L)...
    +4*(L-0.5)*(J-L+0.5)*(J+L+0.5))^2/...
    (2*J*(J+1)*CP(J,Y1,L)*CP(J,Y2,L));

end

% computation of the energy for OH A2Sigma to OH X2Pi
% constant values from Goldman Gilis (1980)

function [E1, E2, E3] = OH_Energy_00(J, p)

x = J+0.5;
y = ((J-0.5)*(J+1.5))^0.5;

H = zeros(6,6);  

B_s = 16.9258978;
D_s = 2.0396682E-03;
H_s = 97.7612E-09;
gamma_s = -7.8395E-03;
nu_o = 32402.056230;
A = -139.228325;
B_p = 18.5497354;
D_p = 1.907852E-03;
H_p = 0.1239836E-06;
A_D = -0.72304E-03;
ALp = -151.9226212;
BLp = 25.0435440;
DLp = 2.6923334E-03;
ADLp = 8.055637E-03;
HLp = 0.166605E-06;


H(1,1) = B_s*x^2+0.25*gamma_s+nu_o-D_s*(x^4+x^2)+H_s*(x^6+3*x^4);
H(1, 2) = B_s*x - 2*D_s*x^3 + H_s*(3*x^5+x^3);
H(2,2) = H(1,1);
H(2,1) = H(1,2);
H(3,3) = B_p*x^2-0.5*A-D_p*(x^4+x^2-1)-A_D*x^2+H_p*(x^6+3*x^4-5*x^2+2);
H(3,1) = BLp+0.5*ALp-DLp*(4*x^2-1)+0.5*ADLp*(2*x^2-1)+HLp*(9*x^4-3*x^2+1);
H(1,3) = H(3,1);
H(4,2) = H(3,1);
H(2,4) = H(4,2);
H(3,2) = BLp*x-DLp*x*(2*x^2+1)+HLp*x*(3*x^4+6*x^2-2);
H(2,3) = H(3,2);
H(4,1) = H(2,3);
H(1,4) = H(4,1);
H(4,4) = H(3,3);
H(5,5) = B_p*(x^2-2)+0.5*A-D_p*(x^4-3*x^2+3)+A_D*(x^2-2)+H_p*(x^6-3*x^4+5*x^2-4);
H(1,5) = BLp*y - DLp*y*(2*x^2-1) + ADLp*y + HLp*y*(3*x^4+1);
H(5,1) = H(1,5);
H(6,2) = H(5,1);
H(2,6) = H(6,2);
H(5,2) = -2*DLp*x*y+3*HLp*x*(2*x^2-1)*y;
H(2,5) = H(5,2);
H(1,6) = H(2,5);
H(6,1) = H(1,6);
H(3,5) = B_p*y-2*D_p*y*(x^2-1)+H_p*y*(3*x^4-5*x^2+3);
H(5,3) = H(3,5);
H(6,4) = H(5,3);
H(4,6) = H(6,4);
H(6,6) = H(5,5);

HK = zeros(3,3);
HK(1,1) = H(1,1)+p*H(1,2);
HK(2,2) = H(3,3);
HK(3,3) = H(5,5);
HK(2,1) = H(1,3)+p*H(2,3);
HK(1,2) = HK(2,1);
HK(1,3) = H(1,5)+p*H(2,5);
HK(3,1) = HK(1,3);
HK(2,3) = H(3,5);
HK(3,2) = HK(2,3);


E = eig(HK);
E3 = E(1);
E2 = E(2);
E1 = E(3);

end
