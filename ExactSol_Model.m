
%%%%%%%%%%%%%%%%%%%%%% Solver for the Element %%%%%%%%%%%%
function [uu, vv, duudx, duudy, dvvdx, dvvdy, Eps11,Eps12, Eps22,Sigma11, Sigma12, Sigma22] = ExactSol_Model(xx, yy, k1, kappa, mu, lambda)

%Get the exact displacement, strain and srresses, according to the
%model-1-crack test case

% Get polar coordinates (rr, th) out of (xx, yy).
rr = sqrt (xx^2 + yy^2);               % crack tips at (0,0)
drrdx = xx/rr;                               % cosin theta  x
drrdy = yy/rr;                               % cosin theta  y

dthdx = -yy/(xx^2 + yy^2);         % dtheta/dx

if xx == 0;
    dthdy = 0;
else
    dthdy = 1/(xx + yy*yy/xx);
end

if xx == 0
    if yy == 0
        disp('Theta can not be determined by (x,y)=(0,0)');
        th = 0;
    elseif yy >0
        th = 0.5* pi;
    else
        th = 1.5*pi;
    end
    
elseif yy == 0
if  xx>0
    th = 0;
else
    th = pi;
end
elseif xx>0
    th = atan(yy/xx);
elseif xx<0
    th = pi + atan(yy/xx);
else
    error ('Internal Error!')
end

if th >pi               % This is important due to the multiplication with 1/2 blow
    th = th - 2* pi;
end

%%%%%%%%%%%%%%%%%%Get Exact Solution%%%%%%%%%%%%

uu = k1/(2*mu)*sqrt(rr/(2*pi))*cos(0.5*th)*(kappa-1+2*sin(0.5*th)*sin(0.5*th));                                 % ?????????
duudrr = 1/8*k1*2^(1/2)*cos(1/2*th)*(kappa+1-2*cos(1/2*th)^2)/mu/pi^(1/2)/rr^(1/2);                 %duudrr
duudth = -1/8*k1*2^(1/2)/pi^(1/2)*rr^(1/2)*sin(1/2*th)*(kappa+1-6*cos(1/2*th)^2)/mu;
duudx = duudrr * drrdx + duudth * dthdx;
duudy = duudrr * drrdy + duudth * dthdy;

vv = k1/(2*mu)*sqrt(rr/(2*pi))*sin(0.5*th)*(kappa+1-2*cos(0.5*th)*cos(0.5*th));
dvvdrr = 1/8*k1/mu*2^(1/2)/pi^(1/2)/rr^(1/2)*sin(1/2*th)*(kappa+1-2*cos(1/2*th)^2);
dvvdth = 1/8*k1*2^(1/2)/pi^(1/2)*rr^(1/2)*cos(1/2*th)*(kappa+5-6*cos(1/2*th)^2)/mu;
dvvdx = dvvdrr * drrdx + dvvdth * dthdx;
dvvdy = dvvdrr * drrdy + dvvdth * dthdy;

Sigma11 = k1/sqrt(2*pi*rr)*cos(0.5*th)*(1-sin(0.5*th)*sin(1.5*th));
Sigma22 = k1/sqrt(2*pi*rr)*cos(0.5*th)*(1+sin(0.5*th)*sin(1.5*th));
Sigma12 = k1/sqrt(2*pi*rr)*cos(0.5*th)*(sin(0.5*th)*cos(1.5*th));

% Sigma = [Sigma11 Sigma12; Sigma12 Sigma22];
% Eps = -0.25*lambda/(mu*(lambda+mu))*trace(Sigma)*[1 0; 0 1] + 1/(2*mu)*Sigma;
% Eps11 = Eps(1,1);
% Eps12 = Eps(1,2);
% Eps22 = Eps(2,2);

Eps11 = duudx;
Eps12 = 0.5 * (duudy + dvvdx);
Eps22 = dvvdy;
end



