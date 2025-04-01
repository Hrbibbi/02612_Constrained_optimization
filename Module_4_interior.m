clear all; close all;
clc
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

function InteriorPoint(H,g,A,b,C,d,x0)
    %---------------------------------
    % Initialization step
    %---------------------------------
    s = ones(size(d)); %maybe check if these are good
    z = ones(size(d)); 
    x = x0;
    y = ones(size(b));
    
    [epsilonL,epsilonA,epsilonC,epsilonmu] = deal(1e-4,1e-4,1e-4,1e-4);
    neq = size(A,2);
    nineq = size(C,2);

    %---------------------------------
    % Iteration step
    %---------------------------------
    for k = 0:101
        %compute residuals
        [rL,rA,rC,rSZ,mu] = computeResidual(H,g,A,b,C,d,x,y,s,z);
        % iteraton print
        fprintf('Iteration %3i, x1 = %7.3f , x2 = %7.3f, |r_L| = %6.4e, |r_A| = %6.4e, |r_C| = %6.4e, |mu| = %6.4e \n', k, x(1), x(2), norm(rL), norm(rA), norm(rC),abs(mu))
        % check if we have converged
        if norm(rL)<epsilonL && norm(rA)<epsilonA && norm(rC)<epsilonC && abs(mu)< epsilonmu
            fprintf('Converged'); break
        end
        if k>=100
            fprintf('max iter reached'); break
        end

        %---------------------------------
        % Precomputation step good
        %---------------------------------
        SinvZ = z./s; ZinvS = s./z;
        % Setup LDL and factorize
        Hbar = H+C*diag(SinvZ)*C';
        [L,D,P] = ldl([Hbar -A;-A',zeros(size(A,2))]);
        %disp(A)
        %---------------------------------
        % Affine Direction step 
        %---------------------------------
        rBarL = rL-C*diag(SinvZ)*(rC-rSZ./z);
        temp = P*(L'\(D\(L\(P'*[-rBarL;-rA]))));
        disp(temp)
        deltaXaff = temp(1:length(x));
        deltaYaff = temp(1+length(x):end);
        deltaZaff = -diag(SinvZ)*C'*deltaXaff+diag(SinvZ)*(rC-rSZ./z);
        deltaSaff = -rSZ./z-ZinvS.*deltaZaff;

        %---------------------------------
        %  First alpha step
        %---------------------------------
        alpha_aff = [-z./deltaZaff;-s./deltaSaff];
        alpha_aff = min(alpha_aff([deltaZaff<0;deltaSaff<0]));

        %---------------------------------
        % Duality Gap
        %---------------------------------
        mu_aff = mean((z+alpha_aff*deltaZaff).*(s+alpha_aff*deltaSaff));
        sigma = (mu_aff/mu)^3;

        %---------------------------------
        % Correction step
        %---------------------------------
        rSZ_bar = (rSZ+deltaSaff.*deltaZaff-sigma*mu);
        rBarL = rL-C*(SinvZ.*(rC-rSZ_bar./z));
        temp = P*(L'\(D\(L\(P'*[-rBarL;-rA])))); %was rA not -rA before
        deltaX = temp(1:length(x));
        deltaY = temp(1+length(x):end);
        deltaZ = -diag(SinvZ)*C'*deltaX+SinvZ.*(rC-rSZ_bar./z);
        deltaS = -rSZ_bar./z-ZinvS.*deltaZ;

        %---------------------------------
        % Second alpha step
        %---------------------------------
        alpha = [-z./deltaZ;-s./deltaS];
        alpha = min(alpha([deltaZ<0;deltaS<0]));

        %---------------------------------
        % Update step
        %---------------------------------
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
A = [1;0] %[1;1];
b = [1]   %[1];
C = [1 -2; -1 -2; -1 2; 1 0; 0 1]';
d = -[2;6;2;0;0];
x0 = [2;0];

contourplt(H,g,C',-d,[0,4,0,4],20)


InteriorPoint(H,g,A,b,C,d,x0);