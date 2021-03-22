% artnetwork.m (1-30-2006) ARTPHONE model of word processing

function dy = artphone(t,y,Mask,IN,ipi,ipa,ipb)
% Main Parameter values: 
% alpha and delta = self-inhibition (decay) parameters
% beta and eta = ceiling value for chunk activation level
% mu and lambda = control connection weight dynamics (not linear)
% kappa = top-down exciation rate parameter
% gamma = threshold for nodes to initiate bottom-up excitation
delta = 0.28; alpha = .5; beta = 1; 
mu = 2.2; lambda = 0.10; g = 1.1;
kappa=.20;gamma=.1;eta=0.1;

%frequency settings for phonemes
pi=1; pn=.999; pc=.97; pa=.8; pb=.6; pr=.6; pt=.55;

a=ipa; %slope parameter for biphones
prob=[0.0144 0.0122 0.0095 0.0059 0.0059 0.0050 0.0041 0.0039]'; %default biphone probs
pwt=a*log(prob)+1-a*log(0.0144);
pan=pwt(1,1);pca=pwt(2,1);pin=pwt(3,1);pba=pwt(4,1);pat=pwt(5,1);pra=pwt(6,1);pbi=pwt(7,1);pta=pwt(8,1);

b=ipb; %slope parameter for words
prob=[0.0133 0.0097 0.0092 0.0068]'; %default word node probs
pnwt=b*log(prob)+1-b*log(.0133);

pfreq=[1 .943932 .885289 .885289]'; % rescaled word frequency values for the four words.

% neighborhood and frequency values are combined.
pcan=pnwt(1,1) * pfreq(1,1); pran=pnwt(2,1) * pfreq(2,1) ;ptan=pnwt(3,1) * pfreq(3,1); pbin=pnwt(4,1) * pfreq(4,1);

%--- w+'s ---
y1m = max([(y(1)-gamma) 0]);y2m = max([(y(2)-gamma) 0]);y3m = max([(y(3)-gamma) 0]);y4m = max([(y(4)-gamma) 0]);
y5m = max([(y(5)-gamma) 0]);y6m = max([(y(6)-gamma) 0]);y7m = max([(y(7)-gamma) 0]);

%--- inputs ---
switch ipi
 case 'TAN'
%TAN - dense word
%%{
Ia=((t>2)&(t<6))*0.1;
In=((t>4)&(t<8))*0.1;
Ic=((t>0)&(t<4))*0.0;
Ii=((t>0)&(t<4))*0.0;
Ib=((t>0)&(t<4))*0.0;
It=((t>0)&(t<4))*0.1;
Ir=((t>0)&(t<4))*0.0;
%}
 case 'BIN'
%BIN - sparse word
%%{
Ia=((t>0)&(t<4))*0.0;
In=((t>4)&(t<8))*0.1;
Ic=((t>0)&(t<4))*0.0;
Ii=((t>2)&(t<6))*0.1;
Ib=((t>0)&(t<4))*0.1;
It=((t>0)&(t<4))*0.0;
Ir=((t>0)&(t<4))*0.0;
%}
  case 'BAN'
%BAN - high prob pword
%%{
Ia=((t>2)&(t<6))*0.1;
In=((t>4)&(t<8))*0.1;
Ic=((t>0)&(t<4))*0.0;
Ii=((t>0)&(t<4))*0.0;
Ib=((t>0)&(t<4))*0.1;
It=((t>0)&(t<4))*0.0;
Ir=((t>0)&(t<4))*0.0;
%}
  case 'BAT'
%BAT - lowest prob pword
%%{
Ia=((t>2)&(t<6))*0.1;
In=((t>4)&(t<8))*0.0;
Ic=((t>0)&(t<4))*0.0;
Ii=((t>0)&(t<4))*0.0;
Ib=((t>0)&(t<4))*0.1;
It=((t>4)&(t<8))*0.1;
Ir=((t>0)&(t<4))*0.0;
%}
end;
% ART network -----------
%phoneme INPUTS in working memory --------------
% Equation 1 in Grossberg et al (1997) except that kappa added to "uz" to
% modulate top-down feedback from sublexical (biphone) nodes
dy=[g*((beta-y(1))*(Ia + kappa*(y(8)*y(28)+y(15)*y(42)+y(16)*y(44)+y(18)*y(46)+y(19)*y(48)+y(20)*y(50)+...
        y(22)*y(52)+y(23)*y(74)+y(24)*y(76)+y(25)*y(78))) - y(1)*alpha);% input 'a'
g*((beta-y(2))*(In + kappa*(y(9)*y(30)+y(15)*y(54)+y(17)*y(56)+y(23)*y(80)+y(24)*y(82)+y(25)*y(84)+y(26)*y(86))) - y(2)*alpha);% input 'n'
g*((beta-y(3))*(Ic + kappa*(y(10)*y(32)+y(16)*y(58)+y(23)*y(88))) - y(3)*alpha);% input 'c'
g*((beta-y(4))*(Ii + kappa*(y(11)*y(34)+y(17)*y(60)+y(21)*y(62)+y(26)*y(90))) - y(4)*alpha);% input 'i'
g*((beta-y(5))*(Ib + kappa*(y(12)*y(36)+y(18)*y(64)+y(21)*y(66)+y(26)*y(92))) - y(5)*alpha);% input 'b'
g*((beta-y(6))*(It + kappa*(y(13)*y(38)+y(19)*y(68)+y(22)*y(70)+y(25)*y(94))) - y(6)*alpha);% input 't'
g*((beta-y(7))*(Ir + kappa*(y(14)*y(40)+y(20)*y(72)+y(24)*y(96))) - y(7)*alpha);% input 'r'

% Sublexical PHONEME nodes
% Equation 3t with masking.
g*((beta-y(8))*y1m*pa*y(27) - delta*y(8) - Mask*(y(15)+y(16)+y(18)+y(19)+y(20)+y(22)+y(23)+y(24)+y(25))); % sublexical 'a'
g*((beta-y(9))*y2m*pn*y(29) - delta*y(9) - Mask*(y(15)+y(17)+y(23)+y(24)+y(25)+y(26))); % sublexical 'n'
g*((beta-y(10))*y3m*pc*y(31) - delta*y(10) - Mask*(y(16)+y(23))); % sublexical 'c'
g*((beta-y(11))*y4m*pi*y(33) - delta*y(11) - Mask*(y(17)+y(21)+y(26))); % sublexical 'i'
g*((beta-y(12))*y5m*pb*y(35) - delta*y(12) - Mask*(y(18)+y(21)+y(26))); % sublexical 'b'
g*((beta-y(13))*y6m*pt*y(37) - delta*y(13) - Mask*(y(19)+y(22)+y(25))); % sublexical 't'
g*((beta-y(14))*y7m*pr*y(39) - delta*y(14) - Mask*(y(20)+y(24))); % sublexical 'r'

% Sublexical BIPHONE nodes
% Equation 3, but expanded to include top-down excitation (kappa) and
% inhibition (latteral and masking).
g*((beta-y(15))*(y1m*pan*y(41)+y2m*pan*y(53)) - delta*y(15) - Mask*(y(23)+y(24)+y(25))-IN*(y(17)+y(19))); % sublexical 'an'
g*((beta-y(16))*(y1m*pca*y(43)+y3m*pca*y(57)) - delta*y(16) - Mask*y(23)-IN*(y(18)+y(20)+y(22))); % sublexical 'ca'
g*((beta-y(17))*(y2m*pin*y(55)+y4m*pin*y(59)) - delta*y(17) - Mask*y(26)-IN*y(15)); % sublexical 'in'
g*((beta-y(18))*(y1m*pba*y(45)+y5m*pba*y(63)) - delta*y(18) - IN*(y(16)+y(20)+y(21)+y(22))); % sublexical 'ba'
g*((beta-y(19))*(y1m*pat*y(47)+y6m*pat*y(67)) - delta*y(19) - IN*y(15)); % sublexical 'at'
g*((beta-y(20))*(y1m*pra*y(49)+y7m*pra*y(71)) - delta*y(20) - Mask*y(24)-IN*(y(16)+y(18)+y(22))); % sublexical 'ra'
g*((beta-y(21))*(y4m*pbi*y(61)+y5m*pbi*y(65)) - delta*y(21) - Mask*y(26)-IN*y(18)); % sublexical 'bi'
g*((beta-y(22))*(y1m*pta*y(51)+y6m*pta*y(69)) - delta*y(22) - Mask*y(25)-IN*(y(16)+y(18)+y(20))); % sublexical 'ta'

%Sublexical WORD nodes
% Equation 3, but expanded to include lateral inhibition (IN).
g*((beta-y(23))*(y1m*pcan*y(73)+y2m*pcan*y(79)+y3m*pcan*y(87)) - delta*y(23) - IN*(y(24)+y(25))); % lexical 'can'
g*((beta-y(24))*(y1m*pran*y(75)+y2m*pran*y(81)+y7m*pran*y(95)) - delta*y(24) - IN*(y(23)+y(25))); % lexical 'ran'
g*((beta-y(25))*(y1m*ptan*y(77)+y2m*ptan*y(83)+y6m*ptan*y(93)) - delta*y(25) - IN*(y(23)+y(24))); % lexical 'tan'
g*((beta-y(26))*(y2m*pbin*y(85)+y4m*pbin*y(89)+y5m*pbin*y(91)) - delta*y(26)); % lexical 'bin'

%Resonance connections
% Even-numbered y vectors are bottom-up; odd-numbered y vectors are top-down.
% lambda and mu control connection-weight dynamcs. See equations 4, 5, 6.
% --- input to phoneme connections
eta*(1-y(27))-(lambda*y1m+mu*y1m^2)*y(27);eta*(1-y(28))-(lambda*y(8)+mu*y(8)^2)*y(28); % bottom-up & top-down connections
eta*(1-y(29))-(lambda*y2m+mu*y2m^2)*y(29);eta*(1-y(30))-(lambda*y(9)+mu*y(9)^2)*y(30); % bottom-up & top-down connections
eta*(1-y(31))-(lambda*y3m+mu*y3m^2)*y(31);eta*(1-y(32))-(lambda*y(10)+mu*y(10)^2)*y(32); % bottom-up & top-down connections
eta*(1-y(33))-(lambda*y4m+mu*y4m^2)*y(33);eta*(1-y(34))-(lambda*y(11)+mu*y(11)^2)*y(34); % bottom-up & top-down connections
eta*(1-y(35))-(lambda*y5m+mu*y5m^2)*y(35);eta*(1-y(36))-(lambda*y(12)+mu*y(12)^2)*y(36); % bottom-up & top-down connections
eta*(1-y(37))-(lambda*y6m+mu*y6m^2)*y(37);eta*(1-y(38))-(lambda*y(13)+mu*y(13)^2)*y(38); % bottom-up & top-down connections
eta*(1-y(39))-(lambda*y7m+mu*y7m^2)*y(39);eta*(1-y(40))-(lambda*y(14)+mu*y(14)^2)*y(40); % bottom-up & top-down connections----
% --- input to biphone connections
eta*(1-y(41))-(lambda*y1m+mu*y1m^2)*y(41);eta*(1-y(42))-(lambda*y(15)+mu*y(15)^2)*y(42); % bottom-up & top-down connections
eta*(1-y(43))-(lambda*y1m+mu*y1m^2)*y(43);eta*(1-y(44))-(lambda*y(16)+mu*y(16)^2)*y(44); % bottom-up & top-down connections
eta*(1-y(45))-(lambda*y1m+mu*y1m^2)*y(45);eta*(1-y(46))-(lambda*y(18)+mu*y(18)^2)*y(46); % bottom-up & top-down connections
eta*(1-y(47))-(lambda*y1m+mu*y1m^2)*y(47);eta*(1-y(48))-(lambda*y(19)+mu*y(19)^2)*y(48); % bottom-up & top-down connections
eta*(1-y(49))-(lambda*y1m+mu*y1m^2)*y(49);eta*(1-y(50))-(lambda*y(20)+mu*y(20)^2)*y(50); % bottom-up & top-down connections
eta*(1-y(51))-(lambda*y1m+mu*y1m^2)*y(51);eta*(1-y(52))-(lambda*y(22)+mu*y(22)^2)*y(52); % bottom-up & top-down connections
eta*(1-y(53))-(lambda*y2m+mu*y2m^2)*y(53);eta*(1-y(54))-(lambda*y(15)+mu*y(15)^2)*y(54); % bottom-up & top-down connections
eta*(1-y(55))-(lambda*y2m+mu*y2m^2)*y(55);eta*(1-y(56))-(lambda*y(17)+mu*y(17)^2)*y(56); % bottom-up & top-down connections
eta*(1-y(57))-(lambda*y3m+mu*y3m^2)*y(57);eta*(1-y(58))-(lambda*y(16)+mu*y(16)^2)*y(58); % bottom-up & top-down connections
eta*(1-y(59))-(lambda*y4m+mu*y4m^2)*y(59);eta*(1-y(60))-(lambda*y(17)+mu*y(17)^2)*y(60); % bottom-up & top-down connections
eta*(1-y(61))-(lambda*y4m+mu*y4m^2)*y(61);eta*(1-y(62))-(lambda*y(21)+mu*y(21)^2)*y(62); % bottom-up & top-down connections
eta*(1-y(63))-(lambda*y5m+mu*y5m^2)*y(63);eta*(1-y(64))-(lambda*y(18)+mu*y(18)^2)*y(64); % bottom-up & top-down connections
eta*(1-y(65))-(lambda*y5m+mu*y5m^2)*y(65);eta*(1-y(66))-(lambda*y(21)+mu*y(21)^2)*y(66); % bottom-up & top-down connections
eta*(1-y(67))-(lambda*y6m+mu*y6m^2)*y(67);eta*(1-y(68))-(lambda*y(19)+mu*y(19)^2)*y(68); % bottom-up & top-down connections
eta*(1-y(69))-(lambda*y6m+mu*y6m^2)*y(69);eta*(1-y(70))-(lambda*y(22)+mu*y(22)^2)*y(70); % bottom-up & top-down connections
eta*(1-y(71))-(lambda*y7m+mu*y7m^2)*y(71);eta*(1-y(72))-(lambda*y(20)+mu*y(20)^2)*y(72); % bottom-up & top-down connections----
% --- input to lexical unit connections
eta*(1-y(73))-(lambda*y1m+mu*y1m^2)*y(73);eta*(1-y(74))-(lambda*y(23)+mu*y(23)^2)*y(74); % bottom-up & top-down connections
eta*(1-y(75))-(lambda*y1m+mu*y1m^2)*y(75);eta*(1-y(76))-(lambda*y(24)+mu*y(24)^2)*y(76); % bottom-up & top-down connections
eta*(1-y(77))-(lambda*y1m+mu*y1m^2)*y(77);eta*(1-y(78))-(lambda*y(25)+mu*y(25)^2)*y(78); % bottom-up & top-down connections
eta*(1-y(79))-(lambda*y2m+mu*y2m^2)*y(79);eta*(1-y(80))-(lambda*y(23)+mu*y(23)^2)*y(80); % bottom-up & top-down connections
eta*(1-y(81))-(lambda*y2m+mu*y2m^2)*y(81);eta*(1-y(82))-(lambda*y(24)+mu*y(24)^2)*y(82); % bottom-up & top-down connections
eta*(1-y(83))-(lambda*y2m+mu*y2m^2)*y(83);eta*(1-y(84))-(lambda*y(25)+mu*y(25)^2)*y(84); % bottom-up & top-down connections
eta*(1-y(85))-(lambda*y2m+mu*y2m^2)*y(85);eta*(1-y(86))-(lambda*y(26)+mu*y(26)^2)*y(86); % bottom-up & top-down connections
eta*(1-y(87))-(lambda*y3m+mu*y3m^2)*y(87);eta*(1-y(88))-(lambda*y(23)+mu*y(23)^2)*y(88); % bottom-up & top-down connections
eta*(1-y(89))-(lambda*y4m+mu*y4m^2)*y(89);eta*(1-y(90))-(lambda*y(26)+mu*y(26)^2)*y(90); % bottom-up & top-down connections
eta*(1-y(91))-(lambda*y5m+mu*y5m^2)*y(91);eta*(1-y(92))-(lambda*y(26)+mu*y(26)^2)*y(92); % bottom-up & top-down connections
eta*(1-y(93))-(lambda*y6m+mu*y6m^2)*y(93);eta*(1-y(94))-(lambda*y(25)+mu*y(25)^2)*y(94); % bottom-up & top-down connections
eta*(1-y(95))-(lambda*y7m+mu*y7m^2)*y(95);eta*(1-y(96))-(lambda*y(24)+mu*y(24)^2)*y(96)]; % bottom-up & top-down connections------- 
