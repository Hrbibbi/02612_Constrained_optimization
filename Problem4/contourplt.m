function contourplt(obj,eqcon,ineqcon,cl,cu,xl,xu,lim,npoint);
    [x1,x2] = meshgrid(linspace(lim(1),lim(2),npoint),linspace(lim(3),lim(4),npoint));
    x = [x1(:),x2(:)]';
    f = zeros(size(x,2),1);
    neqcon = length(eqcon([0,0]));
    nineqcon = length(ineqcon([0,0]));

    ceq = zeros(size(x,2),neqcon);
    cineq = zeros(size(x,2),nineqcon);
    box = zeros(size(x,2),1);
    for i =1:size(x,2)
        f(i) = obj(x(:,i));
        ceq(i,:) = eqcon(x(:,i));
        cineq(i,:) = ineqcon(x(:,i));
        box(i) = all([xl<=x(:,i); xu>=x(:,i)]);
    end
    
    % xbound = ones(size(x,2));
    % xbound(x(1)<xl(1)) =
    

    figure(1)
    clf
    %% Plot objective function
    fval = reshape(f,size(x1));
    contour(x1,x2,fval,60)
    colormap('turbo')
    colorbar
    hold on;
    %% Plot equality constraints
    for i = 1:neqcon
        ceqval = reshape(ceq(:,i),size(x1));
        contour(x1,x2,ceqval,[-0.01 0.01],'red')
    end
    %% Plot inequality constraints
    for i = 1:nineqcon
        cineqval_lower = reshape(cineq(:,i),size(x1))-cl;
        cineqval_upper = reshape(cineq(:,i),size(x1))-cu;
        [c,h] = contourf(x1,x2,-cineqval_lower,[0 0.001],'black');
        h.FaceAlpha = 0.2;  
        h.FaceColor = 'black';
        [c,h] = contourf(x1,x2,cineqval_upper,[0 0.001],'black');
        h.FaceAlpha = 0.2;  
        h.FaceColor = 'black';
    end
    %% Plot box constraints
    boxval = reshape(box,size(x1));
    [c,h] = contourf(x1,x2,-boxval,[-0.001,0.001],'black');
    h.FaceAlpha = 0.2; 
    h.FaceColor = 'black';
    

    xlim([lim(1:2)])
    ylim([lim(3:4)])
end
