% psp.m (1-30-2006) PSP analysis of ART network
% Creates the PSP plots of inhibition by masking.
% This program is an expanded version of Fourway.m

Na=60;
Nb=60;
a = .3; 
b = .017; 
pmask = (0.005:.005:.3)';
pin = (0.0025:.0025:.15)';
zz = zeros(60,60);
yini = [zeros(26,1); ones(70,1)];
tspan = 0:0.2:20;

disp('Total Iterations = 60');

for i=1:60
    disp(num2str(i,3)); 
    for j = 1:60
        inval = 0;
        Mask = pmask(i,1);
        IN = pin(j,1);
        [t,y] = ode45(@artphone,tspan,yini,[],Mask,IN,'TAN',a,b);
        pk1 = y(:,25);
        [a1,b1] = max(pk1);
        
        [t,y] = ode45(@artphone,tspan,yini,[],Mask,IN,'BIN',a,b);
        pk2 = y(:,26);
        [a2,b2] = max(pk2);
        
        [t,y] = ode45(@artphone,tspan,yini,[],Mask,IN,'TAN',a,b);
        pk3 = y(:,15);
        [a3,b3] = max(pk3); %BAT
        
        [t,y] = ode45(@artphone,tspan,yini,[],Mask,IN,'BIN',a,b);
        pk4 = y(:,17);
        [a4,b4] = max(pk4); %BAN
    
        if a1(1,1) > (a2(1,1)+0.02) 
            qq1 = 2;
        elseif abs(a1(1,1)-a2(1,1)) <= 0.02 
            qq1 = 1;
        else
            qq1 = 0;
        end
    
        if a3(1,1) > (a4(1,1)+0.02) 
            qq2 = 2;
        elseif abs(a3(1,1)-a4(1,1)) <= 0.02 
            qq2 = 1;
        else
            qq2 = 0;
        end
   
        if ((a1(1,1) < 0.2) || (a2(1,1) < 0.2) || (a3(1,1) < 0.1) || (a4(1,1) < 0.1)) 
            inval = 1; 
        end
  
        dpattern = 3*qq1+qq2+1;
   
        if (inval == 1) 
            dpattern = 10; 
        end

        zz(i,j) = dpattern/10;
    
    end
end


