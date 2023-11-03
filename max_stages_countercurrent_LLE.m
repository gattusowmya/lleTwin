function stages= max_stages_countercurrent_LLE(F,S)
    
    % with the assumption that the nth raffinate should noth have more than
    % 1.5% of solute, we determine the stages obtained as the maximum
    % stages and then is used as a safety check in case the user has
    % inputted a value too high.
    
    %%
    %F=2000;
    %S=2500;
    xbf = 0.65; % write the value here;
    xcf = 0.35;% write the value here;
    ybs = 0.97;% write the value here;
    ycs = 0.03;
    xcm=0.015; %solute at nth raffinate
    
    %% LLE curve drawing
    
    %B = [0.04 0.05 0.07 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.993 1];
    %C = [0 0.11 0.26 0.375 0.474 0.487 0.468 0.423 0.356 0.274 0.185 0.09 0.001 0];
    
    %data
    % raffinate
    rB=[0.0155 0.017 0.025 0.038 0.06 0.122 0.225 0.354];
    rC=[0 0.0285 0.117 0.205 0.262 0.328 0.346 0.336];
    % extract
    eB=[0.9788 0.9533 0.857 0.735 0.609 0.472 0.354];
    eC=[0 0.0187 0.089 0.173 0.246 0.308 0.336];
    
    
 
    
    
    %% tie line data
    

    tiexc=[0.0285 0.117 0.205 0.262 0.328];
    tieyc=[0.0187 0.089 0.173 0.246 0.308];
%     tiexc = [0.1 0.245 0.426]; % raffinate of solute
%     tiexb = [0.048 0.065 0.133];
%     tieyc = [0.098 0.242 0.409]; % extract of solute
%     tieyb = [0.891 0.736 0.523];
    lt = length(tiexc);
    tie_slope = zeros(1,lt);
    
    
    
    %% interpolation of the points in LLE and drawing the tie lines
    
    B=[rB flip(eB(1:length(eB)-1),2)];
    C=[rC flip(eC(1:length(eC)-1),2)];
    
    B1 = zeros(1,10001);
    C1 = zeros(1,10001);
    for i = 1:10001
    B1(i)=(i-1)/10000;
    end
    
     C1=interp1(B,C,B1);
%     C1=C1(400:end);
%     B1=B1(400:end);

    for i = 1:length(tiexc)
        index1 = find(abs(C1-tiexc(i)) <= 0.0005,1,"first");
        index2 = find(abs(C1-tieyc(i)) <= 0.0005,1,"last");
        XX = [B1(index1) B1(index2)] ;
        YY = [tiexc(i) tieyc(i)];
        
        tie_slope(i) = (YY(2)-YY(1))/(XX(2)-XX(1));
    end
    
    
    %% finding M and Rn
    M = F+S; % overall mass balance
    My = (F*xcf + S*ycs)/(F+S); % solute mass balance
    Mx = (S*ybs)/(F+S);

    
    index3 = find(abs(C1-xcm) <= 0.0005,1,"first");
    xbm=B1(index3);

    
    %% locating the intersection of RnM and FS
    
    
    p1 = polyfit(eB,eC,4);
    p2 = polyfit(rB,rC,4);
    
    
    slopeRnM =(My-xcm)/(Mx-xbm);
    syms x y
    [E1x, E1y] = vpasolve([y==slopeRnM*(x-Mx)+My,y==poly2sym(p1,x)],[x,y] ,[0.3 1; 0.05 0.36]);

    
    %% Finding the intersection of RnS and FE1
    slope1 = (E1y - xcf)/(E1x-0);
    %[ line FE1]
    slope2 = (xcm - 0)/(xbm-ybs);
    % [ line RnS]
    syms x y
    [delx, dely] = vpasolve([y==slope1*(x-E1x)+E1y,y==slope2*(x-xbm)+xcm],[x,y]);

    
    %% finding R1 through tie line interpolation
   
    slope=0;
        if ((E1y>0) & (E1y <= 0.04))
            slope = 0 + (E1y - 0)*tie_slope(1)/(0.04);
        elseif((E1y>0.04) & (E1y <= 0.083))
            slope = tie_slope(1) + (E1y - 0.04)*tie_slope(2)/(0.083 - 0.04);
        elseif ((E1y>0.083) & (E1y <= 0.13))
            slope = tie_slope(2) + (E1y - 0.083)*tie_slope(3)/(0.13 - 0.083);
        elseif ((E1y>0.13) & (E1y <= 0.215))
            slope = tie_slope(3) + (E1y - 0.13)*tie_slope(4)/(0.215 - 0.13);
        elseif ((E1y>0.215 ) & (E1y <= 0.395))
            slope = tie_slope(4) + (E1y - 0.215)*tie_slope(5)/(0.395 - 0.215);
        elseif((E1y > 0.395))
            slope = tie_slope(5) + (E1y - 0.395)*(-0.155)/(0.4 - 0.395);
        end
    
    
    %%%%%%%%%%%%%%%%%%%%% finding the raffinate(R1) point %%%%%%%%%%%%%%%%%%%%
    [xbr,xcr] = vpasolve([y==slope*(x-E1x)+E1y,y==poly2sym(p2,x)],[x,y],[0 0.3;0 0.305]);
    

    
    %% finding the next extract and raffinate compositions
    
    stages=100;
    slope_1=ones(1,stages);
    slope_R=ones(1,stages);
    Ex=ones(1,stages+1);
    Ey=ones(1,stages+1);
    Rx=ones(1,stages+1);
    Ry=ones(1,stages+1);
    
    Ex(1)=E1x;
    Ey(1)=E1y;
    Rx(1)=xbr;
    Ry(1)=xcr;
    ex=E1x;
    ey=E1y;
    i=1;
    while(xcr>xcm)
     % intersection of the raffinate point of previous stage to delta with the equilibrium curve giving Extract point of current stage %%%%%
        slope_R(i)=(Ry(i)-dely)/(Rx(i)-delx);
        rx=Rx(i);
        ry=Ry(i);
    
        [ex,ey] = vpasolve([y==slope_R(i)*(x-rx)+ry,y==poly2sym(p1,x)],[x,y],[0.6 1;0 0.305]);
    
    
      % finding the raffinate %
    
        if ((0 <ey) && (ey <= 0.04))
            slope_1(i) = 0 + (ey - 0)*tie_slope(1)/(0.04);
        elseif((0.04< ey) && (ey <= 0.083))
            slope_1(i) = tie_slope(1) + (ey - 0.04)*tie_slope(2)/(0.083 - 0.04);
        elseif ((0.083 < ey) && (ey <= 0.13))
            slope_1(i) = tie_slope(2) + (ey - 0.083)*tie_slope(3)/(0.13 - 0.083);
        elseif ((0.13 < ey) && (ey <= 0.215))
            slope_1(i) = tie_slope(3) + (ey - 0.13)*tie_slope(4)/(0.215 - 0.13);
        elseif ((0.215 < ey) && (ey <= 0.395))
            slope_1(i) = tie_slope(4) + (ey - 0.215)*tie_slope(5)/(0.395 - 0.215);
        elseif((ey > 0.395))
            slope_1(i) = tie_slope(5) + (ey - 0.395)*(-0.155)/(0.4 - 0.395);
        end
    
        [xbr,xcr] = vpasolve([y==slope_1(i)*(x-ex)+ey,y==poly2sym(p2,x)],[x,y],[0 0.2;0 0.305]);
    
    
    %     Ex(i+1)=ex; Ey(i+1)=ey;
    %     Rx(i+1)=xbr; Ry(i+1)=xcr;
    
        if(xcr>xcm)
        
        Rx(i+1)=xbr; Ry(i+1)=xcr;
        Ex(i+1)=ex; Ey(i+1)=ey;

        i=i+1;
        else
            break;
        end
    
    end
    stages=i-1;
    
    
end
