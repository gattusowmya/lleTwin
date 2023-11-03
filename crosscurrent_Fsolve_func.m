function [xc, yc, percentage_removed, E, R]= crosscurrent_Fsolve_func(S,F,stages) 
    %------------------------------------------------------ calculations
    % function name might be misleading, using fzero instead of fsolve for
    % both countercurrent and crosscurrent
    %stages = 9;
    S = S*ones(stages,1);% write the value here;
    %F = 1000;% write the value here;
    xbf = 0.65; % write the value here;
    xcf = 0.35;% write the value here;
    ybs = 0.97;% write the value here;
    ycs = 0.03; % write the value here;
    % mention co-ordinates of two tie lines
    Mc_initial=xcf*F; % total amount of solute in feed
    Mc_solvent= sum(ycs.*S); % total amount of solute in solvent
    
    
    %% LLE curve drawing
    
    %data
    % raffinate
    rB=[0.0155 0.017 0.025 0.038 0.06 0.122 0.225 0.354];
    rC=[0 0.0285 0.117 0.205 0.262 0.328 0.346 0.336];
    % extract
    eB=[0.9788 0.9533 0.857 0.735 0.609 0.472 0.354];
    eC=[0 0.0187 0.089 0.173 0.246 0.308 0.336];
    

   
    
    
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

    
    % finding the indices corresponding to the given tie lines
    
    tiexb = zeros(1,length(tiexc));
    tieyb = zeros(1,length(tiexc));
    tie_slope = zeros(1,length(tiexc));
    
    %figure(1);
    for i = 1:length(tiexc)
        index1 = find(abs(C1-tiexc(i))< 0.001,1,'first');
        index2 = find(abs(C1-tieyc(i))< 0.001,1,'last');
        XX = [B1(index1) B1(index2)] ;
        YY = [tiexc(i) tieyc(i)];
    
        %%% plotting the tie lines from the data %%%
%         plot([ XX(1) XX(2)],[ YY(1) YY(2)],'c',"Linewidth",1.5)
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
    
    %figure(2);
    for i = 1:stages
    
        M(i) = F+S(i); % overall mass balance
        My(i) = (F*xcf + S(i)*ycs)/(F+S(i)); % solute mass balance
%         p=plot(0:0.01:1,My(i)*ones(101,1));  % a straight line corresponding to y=My(i)
        Mx(i) = (S(i)*ybs)/(F+S(i)); % we can find the x-coordinate through the intersection of y=My(i) and line FS
        
        
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
    
        %Fe = @(x) [y-slope(i)*(x-Mx(i))+My(i);y-poly2sym(p1,x)];
        
    
        r0=[0 0.2];
        xr = fzero(@(x)raffinate_function(x,slope(i),Mx(i),My(i),p2), r0);
        xbr=xr;xcr=slope(i)*(xr-Mx(i))+My(i);
        xc(i)=xcr;
    
        e0=[0.5 1];
        ye = fzero(@(x)extract_function(x,slope(i),Mx(i),My(i),p1), e0);
        ybe=ye;yce=slope(i)*(ye-Mx(i))+My(i);
        yc(i)=yce;
        
        
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
    Mc_solvent= sum(ycs.*S);
    Mc_removed=Mc_extracted-Mc_solvent;
    
    percentage_removed=(Mc_removed/Mc_initial)*100;
    percentage= percentage_removed;
    
    
    for i=1:stages
        extracted=sum(E(1:i).*yc(1:i));
        solvent_quantity= sum(ycs*S(1:i));
        Mc_removed=extracted-solvent_quantity;
    
        percentage_removed(i)=(Mc_removed/Mc_initial)*100;
    end
    

end