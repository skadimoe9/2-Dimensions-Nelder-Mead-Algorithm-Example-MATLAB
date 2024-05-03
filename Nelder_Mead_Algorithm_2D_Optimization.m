%% Nelder-Mead / Simplex Method - Tugas #5 Komputasi Numerik 
% Fadhli Ammar Taqiyuddin Hakim - 2206817396 

clear, clc

STP = 0; 
tolerance = 1e-6; % if needed

% Coefficient
alpha = 1; % For Reflection
beta = 2; % For Expansion 
gamma = 0.5; % For Outside Contraction (Inside Contraction use gamma*(-1)) 
n = 2; % number of dimensional space
iter = 0; % number of iteration

syms f(x1,x2) x1 x2

% Objective Function
x = [x1;x2];

f(x) = @(x1, x2) 12*x1^(2)+4*x2^(2)-12*x1*x2+2*x1;
funct = f(x1, x2); 

% Optimization Variable 
xtrue = [-(1/3);-(1/2)]; 
x0 = [6;12]; % 1st best point
xg1 = [7;11.5]; % 2nd best point
xg2 = [6.7374;13.1515]; % worst point

xbest = x0; % best point for error calculation
xgood = xg1;
xworst = xg2; 

fx0 = f(x0(1),x0(2)); 
fxg1 = f(xg1(1),xg1(2)); 
fxg2 = f(xg2(1),xg2(2)); 

fprintf("Fadhli Ammar Taqiyuddin Hakim - 2206817396\n")
fprintf("Tugas 5 - Komputasi Numerik \n\n")
%% Plot Kontur ============================================================

xa = linspace(-6,7);
xb = linspace(-7,14); 
[X1, X2] = meshgrid(xa,xb); 
fg = f(X1,X2); 

figure 
[Cc,hc] = contourf(X1,X2,fg,... 
    (-50:10:130)); 
clabel(Cc,hc); 
colorbar 
grid 
set(gca,'xtick',[-6:7]) 
set(gca,'ytick',[-7:14]) 
grid 
xlabel(' {\itx}_1 values','FontName','times','FontSize',14); % x-axes label 
ylabel(' {\itx}_2 values','FontName','times','FontSize',14); 
title({'Nelder-Mead - Tugas 5 - Fadhli Ammar T.H.'}, 'FontName', 'times', 'FontSize',14)
grid 
hold on 

bestpoint0 = plot3(x0(1),x0(2),fx0,'bo','LineWidth',2,...
                            'MarkerEdgeColor','k',...
                            'MarkerFaceColor','g',... 
                            'MarkerSize',10); 
goodpoint0 = plot3(xg1(1),xg1(2),fxg1,'bo','LineWidth',2,...
                            'MarkerEdgeColor','k',...
                            'MarkerFaceColor','y',... 
                            'MarkerSize',10); 
worstpoint0 = plot3(xg2(1),xg2(2),fxg2,'bo','LineWidth',2,...
                            'MarkerEdgeColor','k',...
                            'MarkerFaceColor','r',... 
                            'MarkerSize',10); 

    fprintf("xbest = [%f %f], xgood = [%f %f], xworst = [%f %f]\n", x0, xg1, xg2) 
    fprintf("fxbest = %f, fxgood = %f, fxworst = %f\n", fx0, fxg1, fxg2) 
    fprintf("iteration = %d, eabs = -, eapprox = -\n\n", iter)
% End of Plot 
%% Nelder-Mead Start 

while(~STP) 
    xc = (xbest+xgood)/n; % center of the bottom part of the shape 
    xbest_old = xbest;
    xref = xc + alpha*(xc-xworst);  % Reflection
    %fxref = double(f(xref(1),xref(2)));
    
    % Convergence Test 
    while(true)
        % Better than best point? try expand
        if f(xref(1),xref(2)) < f(xbest(1),xbest(2)) 
            xexp = xc + beta*(xc-xworst); 
            % Is the expansion better than the normal reflection? 
            if f(xexp(1),xexp(2)) < f(xref(1),xref(2))  % Yes
                xworst = xgood; 
                xgood = xbest; 
                xbest = xexp; 
                break
            end
            if f(xexp(1),xexp(2)) > f(xref(1),xref(2))  % No 
                xworst = xgood; 
                xgood = xbest; 
                xbest = xref; 
                break
            end 
        end 
    
        % Better than good point? take the reflection point 
        if (f(xgood(1),xgood(2)) > f(xref(1),xref(2))) && (f(xref(1),xref(2)) > f(xbest(1),xbest(2)))
           xworst = xgood; 
           xgood = xref; 
           break
        end 
        
        % Worst than both best and good? Try doing contraction
        if f(xref(1),xref(2)) >= f(xgood(1),xgood(2)) 
            xcon = xc + gamma*(xc-xworst); % Outer Contraction 
            % Good? Take the xcon as the good point
            if f(xcon(1),xcon(2)) < f(xgood(1),xgood(2)) 
               xworst = xgood; 
               xgood = xcon; 
               break
            end 
            % Still not good? Try inner contraction
            if f(xcon(1),xcon(2)) >= f(xgood(1),xgood(2)) 
                xcon = xc - gamma*(xc-xworst); % Inner Contraction
                % Good? Take xcon
                if f(xcon(1),xcon(2)) < f(xgood(1),xgood(2)) 
                   xworst = xgood; 
                   xgood = xcon; 
                   break
                end 
                % No good? Try to shrink the triangle 
                if f(xcon(1),xcon(2)) >= f(xgood(1),xgood(2)) 
                    xworst = xbest + gamma*(xworst - xbest); 
                    xgood = xbest + gamma*(xgood - xbest); 
                    break
                end 
            end 
        end
    end

    iter = iter + 1; 
    eabs = norm(xtrue - xbest); 
    eapprox = norm(xbest_old - xbest);
    fxbest = f(xbest(1),xbest(2)); 
    fxgood = f(xgood(1),xgood(2)); 
    fxworst = f(xworst(1),xworst(2)); 

    if eapprox < tolerance % If only 4 iterations needed, change the if req. to "iter == 4"
        STP = 1; 
    end 
    
    if iter == 1 
        pause 
        delete(bestpoint0), delete(goodpoint0), delete(worstpoint0)
    end

    fprintf("xbest = [%f %f], xgood = [%f %f], xworst = [%f %f]\n", xbest, xgood, xworst) 
    fprintf("fxbest = %f, fxgood = %f, fxworst = %f\n", fxbest, fxgood, fxworst) 
    fprintf("iteration = %d, eabs = %f, eapprox = %f\n\n", iter, eabs, eapprox)
    
    bestpoint = plot3(xbest(1),xbest(2),fxbest,'bo','LineWidth',2,...
                            'MarkerEdgeColor','k',...
                            'MarkerFaceColor','g',... 
                            'MarkerSize',10); 
    goodpoint = plot3(xgood(1),xgood(2),fxgood,'bo','LineWidth',2,...
                                'MarkerEdgeColor','k',...
                                'MarkerFaceColor','y',... 
                                'MarkerSize',10); 
    worstpoint = plot3(xworst(1),xworst(2),fxworst,'bo','LineWidth',2,...
                                'MarkerEdgeColor','k',...
                                'MarkerFaceColor','r',... 
                                'MarkerSize',10); 
    pause
    if STP == 0 
        delete(bestpoint), delete(goodpoint), delete(worstpoint) 
    end 

end