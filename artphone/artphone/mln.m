function dy = mln(t,y,PG,IN,ipi,ipa,ipb)
% Main Parameter values:
% alpha and delta = self-inhibition (decay) parameters
% beta and eta = ceiling value for chunk activation level
% mu and lambda = control connection weight dynamics (not linear)
% kappa = top-down exciation rate parameter
% gamma = threshold for nodes to initiate bottom-up excitation
delta = 0.28; 
alpha = 0.5; 
beta = 1;
mu = 2.2; 
lambda = 0.1; 
g = 1.1; 
kappa = 0.2; 
gamma = PG; 
eta = 0.1;
gammaBi = 0.1;

a=ipa; %slope parameter for biphones 
prob = [0.0144 0.0122 0.0095 0.0059 0.0059 0.0050 0.0041 0.0039]';
pwt = a*log(prob)+1-a*log(0.0144);
pan = pwt(1,1);
pca = pwt(2,1);
pin = pwt(3,1);
pba = pwt(4,1);
pat = pwt(5,1);
pra = pwt(6,1);
pbi = pwt(7,1);
pta = pwt(8,1);

b = ipb; %slope parameter for words
prob = [0.0133 0.0097 0.0092 0.0068]';
pnwt = b*log(prob)+1-b*log(0.01772);
pcan = pnwt(1,1);
pran = pnwt(2,1);
ptan = pnwt(3,1);
pbin = pnwt(4,1);

switch ipi
    case 'TAN'
        Ia = ((t > 2) & (t < 6))*0.1;
        In = ((t > 4) & (t < 8))*0.1;
        Ic = ((t > 0) & (t < 4))*0.0;
        Ii = ((t > 0) & (t < 4))*0.0;
        Ib = ((t > 0) & (t < 4))*0.0;
        It = ((t > 0) & (t < 4))*0.1;
        Ir = ((t > 0) & (t < 4))*0.0;
    case 'BIN'
        Ia = ((t > 2) & (t < 6))*0.0;
        In = ((t > 4) & (t < 8))*0.1;
        Ic = ((t > 0) & (t < 4))*0.0;
        Ii = ((t > 2) & (t < 6))*0.1;
        Ib = ((t > 0) & (t < 4))*0.1;
        It = ((t > 0) & (t < 4))*0.0;
        Ir = ((t > 0) & (t < 4))*0.0;
    case 'BAN'  
        Ia = ((t > 2) & (t < 6))*0.1;
        In = ((t > 4) & (t < 8))*0.1;
        Ic = ((t > 0) & (t < 4))*0.0;
        Ii = ((t > 0) & (t < 4))*0.0;
        Ib = ((t > 0) & (t < 4))*0.1;
        It = ((t > 0) & (t < 4))*0.0;
        Ir = ((t > 0) & (t < 4))*0.0;
    case 'BAT'   
        Ia = ((t > 2) & (t < 6))*0.1;
        In = ((t > 0) & (t < 4))*0.0;
        Ic = ((t > 0) & (t < 4))*0.0;
        Ii = ((t > 0) & (t < 4))*0.0;
        Ib = ((t > 0) & (t < 4))*0.1;
        It = ((t > 4) & (t < 8))*0.1;
        Ir = ((t > 0) & (t < 4))*0.0;
end

%Set the thresholds at which the phoneme-input and sublexical
%nodes will begin bottom-up excitation.
% See text below equation 3.
y1m = max([(y(1)-gamma) 0]);
y2m = max([(y(2)-gamma) 0]);
y3m = max([(y(3)-gamma) 0]);
y4m = max([(y(4)-gamma) 0]);
y5m = max([(y(5)-gamma) 0]);
y6m = max([(y(6)-gamma) 0]);
y7m = max([(y(7)-gamma) 0]);
% w for sublexical biphones
y8m = max([(y(8)-gammaBi) 0]);
y9m = max([(y(9)-gammaBi) 0]);
y10m = max([(y(10)-gammaBi) 0]);
y11m = max([(y(11)-gammaBi) 0]);
y12m = max([(y(12)-gammaBi) 0]);
y13m = max([(y(13)-gammaBi) 0]);
y14m = max([(y(14)-gammaBi) 0]);
y15m = max([(y(15)-gammaBi) 0]);

%%%%%%%%%%%%%%%%%%%
%%% MLN NETWORK %%%
%%%%%%%%%%%%%%%%%%%
%phoneme INPUTS in working memory --------------
% Equation 1 in Grossberg et al (1997) except that kappa added to "uz" to
% modulate top-down feedback from sublexical (biphone) nodes
dy = [...
    g*((beta-y(1))*(Ia+kappa*(y(8)*y(21)+y(9)*y(25)+y(11)*y(33)+y(12)*y(37)+y(13)*y(41)+y(15)*y(49)))-y(1)*alpha);% input 'a'
    g*((beta-y(2))*(In+kappa*(y(8)*y(23)+y(10)*y(29)))-y(2)*alpha);% input 'n'
    g*((beta-y(3))*(Ic+kappa*(y(9)*y(27)))-y(3)*alpha);% input 'c'
    g*((beta-y(4))*(Ii+kappa*(y(10)*y(31)+y(14)*y(45)))-y(4)*alpha);% input 'i'
    g*((beta-y(5))*(Ib+kappa*(y(11)*y(35)+y(14)*y(47)))-y(5)*alpha);% input 'b'
    g*((beta-y(6))*(It+kappa*(y(12)*y(39)+y(15)*y(51)))-y(6)*alpha);% input 't'
    g*((beta-y(7))*(Ir+kappa*(y(13)*y(43)))-y(7)*alpha);% input 'r'

% Sublexical BIPHONE nodes
% Equation 3, but expanded to include top-down excitation (kappa) and
% latteral inhibition (IN).
    g*((beta-y(8))*(y1m*pan*y(20)+y2m*pan*y(22))-delta*y(8)+kappa*(y(16)*y(53)+y(17)*y(57)+y(18)*y(61))-IN*(y(10)+y(12))); % sublexical 'an'
    g*((beta-y(9))*(y1m*pca*y(24)+y3m*pca*y(26))-delta*y(9)+kappa*(y(16)*y(55))-IN*(y(11)+y(13)+y(15))); % sublexical 'ca'
    g*((beta-y(10))*(y2m*pin*y(28)+y4m*pin*y(30))-delta*y(10)+kappa*(y(19)*y(65))-IN*y(8)); % sublexical 'in'
    g*((beta-y(11))*(y1m*pba*y(32)+y5m*pba*y(34))-delta*y(11)-IN*(y(9)+y(13)+y(14)+y(15))); % sublexical 'ba'
    g*((beta-y(12))*(y1m*pat*y(36)+y6m*pat*y(38))-delta*y(12)-IN*y(8)); % sublexical 'at'
    g*((beta-y(13))*(y1m*pra*y(40)+y7m*pra*y(42))-delta*y(13)+kappa*(y(17)*y(59))-IN*(y(9)+y(11)+y(15))); % sublexical 'ra'
    g*((beta-y(14))*(y4m*pbi*y(44)+y5m*pbi*y(46))-delta*y(14)+kappa*(y(19)*y(67))-IN*y(11)); % sublexical 'bi'
    g*((beta-y(15))*(y1m*pta*y(48)+y6m*pta*y(50))-delta*y(15)+kappa*(y(18)*y(63))-IN*(y(9)+y(11)+y(13))); % sublexical 'ta'

%lexical WORD nodes
% Equation 3, but expanded to include lateral inhibition (IN).
    g*((beta-y(16))*(y8m*pcan*y(52)+y9m*pcan*y(54))-delta*y(16)-IN*(y(17)+y(18))); % lexical 'can'
    g*((beta-y(17))*(y8m*pran*y(56)+y13m*pran*y(58))-delta*y(17)-IN*(y(16)+y(18))); % lexical 'ran'
    g*((beta-y(18))*(y8m*ptan*y(60)+y15m*ptan*y(62))-delta*y(18)-IN*(y(16)+y(17))); % lexical 'tan'
    g*((beta-y(19))*(y10m*pbin*y(64)+y14m*pbin*y(66))- delta*y(19)); % lexical 'bin'

% Resonance connections
% Even-numbered y vectors are bottom-up; odd-numbered y vectors are top-down.
% lambda and mu control connection-weight dynamcs. See equations 4, 5, 6.
% --- phoneme-to-biphone unit connections
    eta*(1-y(20))-(lambda*y1m+mu*y1m^2)*y(20);
    eta*(1-y(21))-(lambda*y(8)+mu*y(8)^2)*y(21); 
    eta*(1-y(22))-(lambda*y2m+mu*y2m^2)*y(22);
    eta*(1-y(23))-(lambda*y(8)+mu*y(8)^2)*y(23); 
    eta*(1-y(24))-(lambda*y1m+mu*y1m^2)*y(24);
    eta*(1-y(25))-(lambda*y(9)+mu*y(9)^2)*y(25); 
    eta*(1-y(26))-(lambda*y3m+mu*y3m^2)*y(26);
    eta*(1-y(27))-(lambda*y(9)+mu*y(9)^2)*y(27); 
    eta*(1-y(28))-(lambda*y2m+mu*y2m^2)*y(28);
    eta*(1-y(29))-(lambda*y(10)+mu*y(10)^2)*y(29); 
    eta*(1-y(30))-(lambda*y4m+mu*y4m^2)*y(30);
    eta*(1-y(31))-(lambda*y(10)+mu*y(10)^2)*y(31); 
    eta*(1-y(32))-(lambda*y1m+mu*y1m^2)*y(32);
    eta*(1-y(33))-(lambda*y(11)+mu*y(11)^2)*y(33); 
    eta*(1-y(34))-(lambda*y5m+mu*y5m^2)*y(34);
    eta*(1-y(35))-(lambda*y(11)+mu*y(11)^2)*y(35); 
    eta*(1-y(36))-(lambda*y1m+mu*y1m^2)*y(36);
    eta*(1-y(37))-(lambda*y(12)+mu*y(12)^2)*y(37); 
    eta*(1-y(38))-(lambda*y6m+mu*y6m^2)*y(38);
    eta*(1-y(39))-(lambda*y(12)+mu*y(12)^2)*y(39); 
    eta*(1-y(40))-(lambda*y1m+mu*y1m^2)*y(40);
    eta*(1-y(41))-(lambda*y(13)+mu*y(13)^2)*y(41); 
    eta*(1-y(42))-(lambda*y7m+mu*y7m^2)*y(42);
    eta*(1-y(43))-(lambda*y(13)+mu*y(13)^2)*y(43); 
    eta*(1-y(44))-(lambda*y4m+mu*y4m^2)*y(44);
    eta*(1-y(45))-(lambda*y(14)+mu*y(14)^2)*y(45); 
    eta*(1-y(46))-(lambda*y5m+mu*y5m^2)*y(46);
    eta*(1-y(47))-(lambda*y(14)+mu*y(14)^2)*y(47); 
    eta*(1-y(48))-(lambda*y1m+mu*y1m^2)*y(48);
    eta*(1-y(49))-(lambda*y(15)+mu*y(15)^2)*y(49); 
    eta*(1-y(50))-(lambda*y6m+mu*y6m^2)*y(50);
    eta*(1-y(51))-(lambda*y(15)+mu*y(15)^2)*y(51);

% --- biphone-to-lexical unit connections
    eta*(1-y(52))-(lambda*y8m+mu*y8m^2)*y(52);
    eta*(1-y(53))-(lambda*y(16)+mu*y(16)^2)*y(53); 
    eta*(1-y(54))-(lambda*y9m+mu*y9m^2)*y(54);
    eta*(1-y(55))-(lambda*y(16)+mu*y(16)^2)*y(55); 
    eta*(1-y(56))-(lambda*y8m+mu*y8m^2)*y(56);
    eta*(1-y(57))-(lambda*y(17)+mu*y(17)^2)*y(57); 
    eta*(1-y(58))-(lambda*y13m+mu*y13m^2)*y(58);
    eta*(1-y(59))-(lambda*y(17)+mu*y(17)^2)*y(59); 
    eta*(1-y(60))-(lambda*y8m+mu*y8m^2)*y(60);
    eta*(1-y(61))-(lambda*y(18)+mu*y(18)^2)*y(61); 
    eta*(1-y(62))-(lambda*y15m+mu*y15m^2)*y(62);
    eta*(1-y(63))-(lambda*y(18)+mu*y(18)^2)*y(63); 
    eta*(1-y(64))-(lambda*y10m+mu*y10m^2)*y(64);
    eta*(1-y(65))-(lambda*y(19)+mu*y(19)^2)*y(65); 
    eta*(1-y(66))-(lambda*y14m+mu*y14m^2)*y(66);
    eta*(1-y(67))-(lambda*y(19)+mu*y(19)^2)*y(67)]; 


