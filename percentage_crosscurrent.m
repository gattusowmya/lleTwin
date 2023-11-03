function percentage = percentage_crosscurrent(stages, S,F)
    %------------------------------------------------------ calculations
    xbf = 0.65; 
    xcf = 0.35;
    ybs = 0.97;
    ycs = 0.03; 
    
    Mc_initial=xcf*F; % total amount of solute in feed
    Mc_solvent= stages*sum(ycs.*S); % total amount of solute in solvent
    
    
    %% LLE curve drawing
    
    %data
    % raffinate
    rB=[0.0155 0.017 0.025 0.038 0.06 0.122 0.225 0.354];
    rC=[0 0.0285 0.117 0.205 0.262 0.328 0.346 0.336];
    % extract
    eB=[0.9788 0.9533 0.857 0.735 0.609 0.472 0.354];
    eC=[0 0.0187 0.089 0.173 0.246 0.308 0.336];
    
    
    figure(1)
    plot(rB,rC,'bo','DisplayName','raffinate data');grid on;
    hold on;
    plot(rB,rC,'r','linewidth',1,'DisplayName','raffinate data');
    plot([0 1 0 0],[0 0 1 0],'k-.','linewidth',1.25)
    plot(eB,eC,'ro','DisplayName','extract data');grid on;
    plot(eB,eC,'b','linewidth',1,'DisplayName','extract data');
    xlabel('xB xC');ylabel('xCyC');
    
    X = [0 0.98]; % x-coordinates of F & S respectively
    Y = [0.35 0.02]; % y-coordinates of F & S respectively 
    plot(X,Y,"m-","linewidth",1.35)
    text(0,0.35,"F")
    text(0.98,0.02,"- S")
    
    
    %% Tie line drawing for n tie lines where n= length(tiexc)
    
    tiexc=[0.0285 0.117 0.205 0.262 0.328];
    tieyc=[0.0187 0.089 0.173 0.246 0.308];
    
    % interpolation of the points in LLE curve
    B=[rB flip(eB(1:length(eB)-1),2)];
    C=[rC flip(eC(1:length(eC)-1),2)];
    
    
    B1= zeros(1,10001);
    C1= zeros(1,10001);
    
    for i= 1:10001
        B1(i)=(i-1)/10000;
        %C1(i) =(interp1(B,C,B1));
    end
     C1=(interp1(B,C,B1));
    figure(2)
    %plot(B1,C1,'b','linewidth',1.00);grid on;hold on;
    %plot([0 1 0 0],[0 0 1 0],'k-.','linewidth',1.25)
    xlabel('xB');ylabel('xC,yC');title('interpolated data - overall process schematic')
    
    % finding the indices corresponding to the given tie lines
    
    tiexb = zeros(1,length(tiexc));
    tieyb = zeros(1,length(tiexc));
    tie_slope = zeros(1,length(tiexc));
    
    figure(1);
    for i = 1:length(tiexc)
        index1 = find(abs(C1-tiexc(i))< 0.001,1,'first');
        index2 = find(abs(C1-tieyc(i))< 0.001,1,'last');
        XX = [B1(index1) B1(index2)] ;
        YY = [tiexc(i) tieyc(i)];
    
        %%% plotting the tie lines from the data %%%
        %plot([ XX(1) XX(2)],[ YY(1) YY(2)],'c',"Linewidth",1.5)
        tie_slope(i) = (YY(2)-YY(1))/(XX(2)-XX(1));
    end
    
    
    
    
    %% finding co-ordinates of M for each stage and the amount of Extract and Raffinate
    R = ones(1,stages);
    E = ones(1,stages);
    xbr = ones(1,stages);
    xcr = ones(1,stages);
    ybe = ones(1,stages);
    yce = ones(1,stages);
    xb = ones(1,stages);
    xc = ones(1,stages);
    yb = ones(1,stages);
    yc = ones(1,stages);
    M = ones(1,stages);
    Mx = ones(1,stages);
    My = ones(1,stages);
    slope= ones(1,stages);
    
    
    p1 = polyfit(eB,eC,6);
    p2 = polyfit(rB,rC,5);
    
    
    syms x y
    
    figure(2);
    for i = 1:stages
    
        M(i) = F+S; % overall mass balance
        My(i) = (F*xcf + S*ycs)/(F+S); % solute mass balance
        %p=plot(0:0.01:1,My(i)*ones(101,1));  % a straight line corresponding to y=My(i)
        Mx(i) = (S*ybs)/(F+S); % we can find the x-coordinate through the intersection of y=My(i) and line FS
        
        
        if ((0 < My(i)) && (My(i) <= 0.04))
            slope(i) = 0 + (My(i) - 0)*tie_slope(1)/(0.04);
        elseif((0.04< My(i)) && (My(i) <= 0.083))
            slope(i) = tie_slope(1) + (My(i) - 0.04)*tie_slope(2)/(0.083 - 0.04);
        elseif ((0.083 < My(i)) && (My(i) <= 0.13))
            slope(i) = tie_slope(2) + (My(i) - 0.083)*tie_slope(3)/(0.13 - 0.083);
        elseif ((0.13 < My(i)) && (My(i) <= 0.215))
            slope(i) = tie_slope(3) + (My(i) - 0.13)*tie_slope(4)/(0.215 - 0.13);
        elseif ((0.215 < My(i)) && (My(i) <= 0.395))
            slope(i) = tie_slope(4) + (My(i) - 0.215)*tie_slope(5)/(0.395 - 0.215);
        elseif((My(i) > 0.395))
            slope(i) = tie_slope(5) + (My(i) - 0.395)*(-0.155)/(0.4 - 0.395);
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% solving for the extract and raffinate compositions in each stage %%%%%%%%%%%%%%%%%%%%% 
        

    
        [xbr,xcr] = vpasolve([y==slope(i)*(x-Mx(i))+My(i),y==poly2sym(p2,x)],[x,y],[0 0.1;0 0.305]); %% raffinate part
        [ybe,yce] = vpasolve([y==slope(i)*(x-Mx(i))+My(i),y==poly2sym(p1,x)],[x,y],[0.6 1;0 0.305]); %% extract part
    
        %%% plotting the mixture point and the extract and raffinate composition point %%%
        xb(i)=xbr(1); xc(i)=xcr(1); yb(i)=ybe(1); yc(i)= yce(1); % saving the vpa solve solution in an array
        %figure(1);
        %plot([double(xb(i)) double(yb(i))],[double(xc(i)) double(yc(i))],"o-","Color",[0,0.25,1],"linewidth",1.35)
        %hold on
        %plot(double(Mx(i)),double(My(i)),"bo")
        text(double(Mx(i)),double(My(i)+0.05),"M")
        
        %%% setting the feed of the next stage as the raffinate of the current stage, likewise the composition of point F %%%
        R(i)=M(i)*(My(i)-yce)/(xcr-yce);
        E(i)=M(i)-R(i);
        F=R(i);
    
        xcf=xcr;
        xbf=xbr;
    
    
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % removal of solute from the feed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %Mc_tot=xcf*F + sum(ycs.*S); % total amount of solute in feed + solvent
    Mc_extracted=sum(E.*yc);
    Mc_removed=Mc_extracted-Mc_solvent;
    
    percentage_removed=(Mc_removed/Mc_initial)*100;
    percentage= percentage_removed;
end