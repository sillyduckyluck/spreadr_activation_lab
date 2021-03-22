% Fourway.m 
%This program runs the ART network once with the fixed set of parameters.
% It then plots the resonance functions for the four conditions in the V&L paper.
% With the current values, the V&L reversal will be plotted.
% Resonance functions for all nodes (words, biphones, phonemes) can be examined by changing y().
a=.3;     %slope parameter for biphones
b=.017;	  %slope parameter for words
zz=zeros(30,30);yini=[zeros(26,1);ones(70,1)];   %nodes initialized to 0, connections to 1

Mask=0.04;IN=0.13;
[t,y]=ode45(@artnetwork,[0:1/1:20],yini,[],Mask,IN,'TAN',a,b);pk1=y(:,25);[a1,b1]=max(pk1);
[t,y]=ode45(@artnetwork,[0:1/1:20],yini,[],Mask,IN,'BIN',a,b);pk2=y(:,26);[a2,b2]=max(pk2);
[t,y]=ode45(@artnetwork,[0:1/1:20],yini,[],Mask,IN,'TAN',a,b);pk3=y(:,15);[a3,b3]=max(pk3);  %BAN
[t,y]=ode45(@artnetwork,[0:1/1:20],yini,[],Mask,IN,'BIN',a,b);pk4=y(:,17);[a4,b4]=max(pk4);  %BAT
plot(t,pk1,'k-o',t,pk2,'k-s',t,pk3,'r:o',t,pk4,'r:s','LineWidth',2);
%plot(t,pk1,'k-o',t,pk2,'k-s',t,pk3,'k:o',t,pk4,'k:s','LineWidth',2);
ylim([-0.1 .75]);
legend('High prob word (tan)','Low prob word (bin)','High prob sublex word (tan)','Low prob sublex word (bin)');
