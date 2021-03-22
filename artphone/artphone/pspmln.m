
a = 0.025; % Slope parameter for biphones 
b = 0.025; % Slope parameter for words 
pin = (0.0025:.0025:.15)';
%pgam = (0.005:.005:.15)';
PG = 0.01;
zz = zeros(60,1);
yini = [zeros(19,1); ones(48,1)]; 
tspan = 0:0.2:20; 

disp('Total Iterations = 60');
for i = 1:60
    disp(num2str(i,3)); 
    inval = 0;
    IN = pin(i,1);
        
    [t,y] = ode45(@mln,tspan,yini,[],PG,IN,'TAN',a,b);
    pk1 = y(:,18);
    [a1,b1] = max(pk1);
    
    [t,y] = ode45(@mln,tspan,yini,[],PG,IN,'BIN',a,b);
    pk2 = y(:,19);[a2,b2]=max(pk2);
    
    [t,y] = ode45(@mln,tspan,yini,[],PG,IN,'TAN',a,b);
    pk3 = y (:,8);
    [a3,b3] = max(pk3); %BAN
    
    [t,y] = ode45(@mln,tspan,yini,[],PG,IN,'BIN',a,b);
    pk4 = y(:,10);
    [a4,b4] = max(pk4); %12 BAT
    
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
    
    zz(i,1) = dpattern/10;
    
    %disp(num2str(10*zz,0));
end



