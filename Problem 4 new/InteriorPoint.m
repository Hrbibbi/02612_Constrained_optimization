
function [x,y,z] = InteriorPoint(H,g,A,b,C,d,x0)
    %---------------------------------
    % Initialization step
    %---------------------------------
    s = ones(size(d));
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
        %fprintf('Iteration %3i, x1 = %7.3f , x2 = %7.3f, |r_L| = %6.4e, |r_A| = %6.4e, |r_C| = %6.4e, |mu| = %6.4e \n', k, x(1), x(2), norm(rL), norm(rA), norm(rC),abs(mu))
        
        % check if we have converged
        if norm(rL)<epsilonL && norm(rA)<epsilonA && norm(rC)<epsilonC && abs(mu)< epsilonmu
            break
        end
        if k>=100
            error('max iter reached');
        end

        %---------------------------------
        % Precomputation step good
        %---------------------------------

        SinvZ = z./s; ZinvS = s./z;
        % Setup LDL and factorize
        Hbar = H+C*diag(SinvZ)*C';
        [L,D,P] = ldl([Hbar -A;-A',zeros(size(A,2))]); %good

        %---------------------------------
        % Affine Direction step 
        %---------------------------------

        %Proposed change to improve speed
        %rBarL = rL- C * ( SinvZ.*(rC-rSZ./z) )
        rBarL = rL-C*diag(SinvZ)*(rC-rSZ./z);

        temp = P*(L'\(D\(L\(P'*[-rBarL;-rA])))); %good
        deltaXaff = temp(1:length(x));
        deltaYaff = temp(1+length(x):end);

        %Proposed change to improve speed
        %deltaZaff = -SinvZ .* C'*deltaXaff + SinvZ .* (rC-rSZ./z)
        deltaZaff = -diag(SinvZ)*C'*deltaXaff+diag(SinvZ)*(rC-rSZ./z);
     
        deltaSaff = -rSZ./z-ZinvS.*deltaZaff;

        %---------------------------------
        %  First alpha step
        %---------------------------------

        alpha_aff = [-z./deltaZaff;-s./deltaSaff]; %good
        num_zeros=alpha_aff([deltaZaff<0;deltaSaff<0]);
        if size(num_zeros,1)==0
            alpha_aff = 1;
        else 
            alpha_aff = min(alpha_aff([deltaZaff<0;deltaSaff<0]));
        end
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
        deltaZ = -SinvZ.*C'*deltaX+SinvZ.*(rC-rSZ_bar./z);
        deltaS = -rSZ_bar./z-ZinvS.*deltaZ;

        %---------------------------------
        % Second alpha step
        %---------------------------------

        alpha = [-z./deltaZ;-s./deltaS];
        num_zeros=alpha([deltaZ<0;deltaS<0]);
        if size(num_zeros,1)==0 
            alpha = 0;
        else
            alpha = min(alpha([deltaZ<0;deltaS<0]));
        end
        %---------------------------------
        % Update step
        %---------------------------------
        alpha_bar = 0.995*alpha;
        %plot([x(1),x(1)+alpha_bar*deltaX(1)],[x(2),x(2)+alpha_bar*deltaX(2)],'-or')
        x = x+alpha_bar*deltaX;
        y = y+alpha_bar*deltaY;
        z = z+alpha_bar*deltaZ;
        s = s+alpha_bar*deltaS;
    end
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
