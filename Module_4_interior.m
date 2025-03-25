clear all; close all;
function contourplt(H,g,A,b,lim,npoint);
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

function [rL,rA,rC,rSZ,mu] = computeResidual(H,g,A,b,C,d,x,y,s,z)
    if isempty(A)
        rA = [];
        rL = H*x+g-C*z;
    elseif isempty(C)
        error("You are and idiot")
    else
        rA = b-A'*x;
        rL = H*x+g-A*y-C*z;
    end
    rC = s+d-C'*x;
    rSZ = s.*z;
    mu = mean(rSZ);
end

function InteriorPoint(H,g,A,b,C,d,x0);
    s = ones(size(d));
    z = ones(size(d)); 
    x = x0;
    y = ones(size(b));
    [epsilonL,epsilonA,epsilonC,epsilonmu] = deal(1e-4,1e-4,1e-4,1e-4);
    neq = size(A,2);
    nineq = size(C,2);
    for k = 0:100
        % Compute residuals
        [rL,rA,rC,rSZ,mu] = computeResidual(H,g,A,b,C,d,x,y,s,z);
        % Check Optimality
        fprintf('Iteration %3i, x1 = %7.3f , x2 = %7.3f, |r_L| = %6.4e, |r_A| = %6.4e, |r_C| = %6.4e, |mu| = %6.4e \n', k, x(1), x(2), norm(rL), norm(rA), norm(rC),abs(mu))
        if norm(rL)<epsilonL && norm(rA)<epsilonA && norm(rC)<epsilonC && abs(mu)< epsilonmu
            fprintf('Converged '); break
        end
        SZ = z./s; ZS = s./z;
        % Setup LDL and factorize
        Hbar = H+C*diag(SZ)*C';
        [L,D,P] = ldl([Hbar -A;-A',zeros(size(A,2))]);
        % Affine Direction
        rBarL = rL-C*diag(SZ)*(rC-rSZ./z);
        temp = P*(L'\(D\(L\(P'*[-rBarL;-rA]))));
        deltaXaff = temp(1:length(rBarL));
        deltaYaff = temp((1+length(rBarL)):(length(rBarL)+length(rA)));
        deltaZaff = -diag(SZ)*C'*deltaXaff+diag(SZ)*(rC-rSZ./z);
        deltaSaff = -rSZ./z-diag(ZS)*deltaZaff;
        alpha_aff = [-z./deltaZaff;-s./deltaSaff];
        alpha_aff = min(alpha_aff([deltaZaff<0;deltaSaff<0]));
        %Duality Gap
        mu_aff = mean((z+alpha_aff*deltaZaff).*(s+alpha_aff*deltaSaff));
        sigma = (mu_aff/mu)^3;
        % Affine Centering Correction direction
        rSZ_bar = (rSZ+deltaSaff.*deltaZaff-sigma*mu);
        rBarL = rL-C*diag(SZ)*(rC-rSZ_bar./z);
        temp = P*(L'\(D\(L\(P'*[-rBarL;rA]))));
        deltaX = temp(1:length(rBarL));
        deltaY = temp((1+length(rBarL)):(length(rBarL)+length(rA)));
        deltaZ = -diag(SZ)*C'*deltaX+diag(SZ)*(rC-rSZ_bar./z);
        deltaS = -rSZ_bar./z-diag(ZS)*deltaZ;
        alpha = [-z./deltaZ;-s./deltaS];
        alpha = min(alpha([deltaZ<0;deltaS<0]));
        alpha_bar = 0.995*alpha;
        plot([x(1),x(1)+alpha_bar*deltaX(1)],[x(2),x(2)+alpha_bar*deltaX(2)],'-or')
        x = x+alpha_bar*deltaX;
        y = y+alpha_bar*deltaY;
        z = z+alpha_bar*deltaZ;
        s = s+alpha_bar*deltaS;   
    end
end

% QP problem
H = 2*eye(2);
g = [-2;-5];
A = [];
b = [];
C = [1 -2; -1 -2; -1 2; 1 0; 0 1]';
d = -[2;6;2;0;0];
x0 = [2;0];

contourplt(H,g,C',-d,[0,4,0,4],20)

InteriorPoint(H,g,A,b,C,d,x0);