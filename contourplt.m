function contourplt(H,g,A,b,lim,npoint)
    [x1,x2] = meshgrid(linspace(lim(1),lim(2),npoint),linspace(lim(3),lim(4),npoint));
    x = [x1(:),x2(:)]';
    for i =1:size(x,2)
        fval(:,i) = g'*x(:,i)+0.5*x(:,i)'*H'*x(:,i);
    end
    fval = reshape(fval,size(x1));
    x1 = reshape(x(1,:),size(x1));
    x2 = reshape(x(2,:),size(x1));
    figure()
    contour(x1,x2,fval)
    hold on;
    for i = 1:length(b)
        if A(i,2) ==0
            x2 = linspace(lim(3),lim(4),npoint);
            x1 = -(A(i,2)*x2+b(i))/A(i,1);
            plot(x1,x2,':k');
        else
            x1 = linspace(lim(1),lim(2),npoint);
            x2 = -(A(i,1)*x1+b(i))/A(i,2);
            plot(x1,x2,':k');
        end
    end
    xlim([lim(1:2)])
    ylim([lim(3:4)])
end