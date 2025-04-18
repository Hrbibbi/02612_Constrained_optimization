function [x] = Interiorpointuppperlower(x0,H,g,C,l,u,dl,du)
    %---------------------------------
    % Initialization step
    %---------------------------------
    
    [yl,yu,zl,zu] = deal(ones(size(l)),ones(size(u)),ones(size(dl)),ones(size(du)));
    [sl,su,tl,tu] = deal(ones(size(l)),ones(size(u)),ones(size(dl)),ones(size(du)));
    x = x0;
    
    [epsilonL,epsilonC,epsilonmu] = deal(1e-4,1e-4,1e-4);
    nineq = size(C,2);

    %---------------------------------
    % Iteration step
    %---------------------------------
    for k = 0:200
        %---------------------------------
        % Residual step
        %---------------------------------
        
        rL=H*x+g+(yu-yl)+C*(zu-zl);
        
        [rc1 ,rc2 ,rc3 ,rc4 ] = deal( sl+l-x, su-u+x, tl+dl-C'*x, tu-du+C'*x );
        rC=[rc1; rc2; rc3; rc4 ];
        
        [rsz1,rsz2,rsz3,rsz4] = deal( sl.*yl, su.*yu, tl.*zl  , tu.*zu );
        
        Z=[yl;yu;zl;zu]; S=[sl;su;tl;tu];

        mu=mean(Z.*S);
        % iteraton print
        fprintf('Iteration %3i, |r_L| = %6.4e, |r_C| = %6.4e, |mu| = %6.4e \n', k, norm(rL), norm(rC),abs(mu))
        % check if we have converged
        if norm(rL)<epsilonL && norm(rC)<epsilonC && abs(mu)< epsilonmu
            fprintf('Converged'); break
        end
        if k>=200
            fprintf('max iter reached'); break
        end

        %---------------------------------
        % Precomputation step
        %---------------------------------
        
        temp1=diag( yl./sl + yu./su );
        temp2=diag( zl./tl + zu./tu );
        Hbar=H+temp1+C*temp2*C';

        R=chol(Hbar);
        
        %---------------------------------
        % Affine Direction step 
        %---------------------------------
        
        y_part=(yu.*rc2-rsz2)./su-(yl.*rc1-rsz1)./sl;
        z_part=(zu.*rc4-rsz4)./tu-(zl.*rc3-rsz3)./tl;
        rL_bar=rL+y_part+C*z_part;
        
        delta_x_affine=R \ ( R' \ (-rL_bar) );
        
        [delta_yl_affine,delta_yu_affine,delta_zl_affine,delta_zu_affine] = deal( ...
           (yl.*(rc1-delta_x_affine)-rsz1)./sl, ...
           (yu.*(rc2+delta_x_affine)-rsz2)./su, ...
           (zl.*(rc3-C'*delta_x_affine)-rsz3)./tl, ...
           (zu.*(rc4+C'*delta_x_affine)-rsz4)./tu  ...
        );
        
        [delta_sl_affine,delta_su_affine,delta_tl_affine,delta_tu_affine] = deal( ...
            -((rsz1+sl.*delta_yl_affine)./yl), ...
            -((rsz2+su.*delta_yu_affine)./yu), ...
            -((rsz3+tl.*delta_zl_affine)./zl), ...
            -((rsz4+tu.*delta_zu_affine)./zu)  ...
        );
        
        deltaZaff=[delta_yl_affine; delta_yu_affine; delta_zl_affine; delta_zu_affine];
        deltaSaff=[delta_sl_affine; delta_su_affine; delta_tl_affine; delta_tu_affine];
        
        %---------------------------------
        %  First alpha step
        %---------------------------------
        
        alpha_aff = [-Z./deltaZaff;-S./deltaSaff];
        num_zeros=alpha_aff([deltaZaff<0;deltaSaff<0]);
        if size(num_zeros,1)==0
            alpha_aff = 1;
        else
            alpha_aff = min(alpha_aff([deltaZaff<0;deltaSaff<0]));
        end
        %---------------------------------
        % Duality Gap
        %---------------------------------
        
        mu_aff = mean((Z+alpha_aff*deltaZaff).*(S+alpha_aff*deltaSaff));
        sigma = (mu_aff/mu)^3;

        %---------------------------------
        % Correction step
        %---------------------------------
        
        [rsz1_bar,rsz2_bar,rsz3_bar,rsz4_bar] = deal(...
          rsz1+delta_sl_affine.*delta_yl_affine-sigma*mu, ...
          rsz2+delta_su_affine.*delta_yu_affine-sigma*mu, ...
          rsz3+delta_tl_affine.*delta_zl_affine-sigma*mu, ...
          rsz4+delta_tu_affine.*delta_zu_affine-sigma*mu ...
        );
        
        y_part=(yu.*rc2-rsz2_bar)./su-(yl.*rc1-rsz1_bar)./sl;
        z_part=(zu.*rc4-rsz4_bar)./tu-(zl.*rc3-rsz3_bar)./tl;
        rL_bar=rL+y_part+C*z_part;

        delta_x= R \ ( R' \ (-rL_bar) );

        [delta_yl,delta_yu,delta_zl,delta_zu] = deal( ...
           (yl.*(rc1-delta_x)-rsz1_bar)./sl, ...
           (yu.*(rc2+delta_x)-rsz2_bar)./su, ...
           (zl.*(rc3-C'*delta_x)-rsz3_bar)./tl, ...
           (zu.*(rc4+C'*delta_x)-rsz4_bar)./tu  ...
        );
        
        [delta_sl,delta_su,delta_tl,delta_tu] = deal( ...
            -((rsz1_bar+sl.*delta_yl)./yl), ...
            -((rsz2_bar+su.*delta_yu)./yu), ...
            -((rsz3_bar+tl.*delta_zl)./zl), ...
            -((rsz4_bar+tu.*delta_zu)./zu)  ...
        );
    
        deltaZ=[delta_yl;delta_yu;delta_zl;delta_zu];
        deltaS=[delta_sl;delta_su;delta_tl;delta_tu];
        %---------------------------------
        % Second alpha step
        %---------------------------------
        
        alpha = [-Z./deltaZ;-S./deltaS];
        num_zeros=alpha([deltaZ<0;deltaS<0]);
        if size(num_zeros,1)==0
            alpha = 1;
        else
            alpha = min(alpha([deltaZ<0;deltaS<0]));
        end
        %---------------------------------
        % Update step
        %---------------------------------
        alpha_bar = 0.995*alpha;
        
        x = x+alpha_bar*delta_x;

        [yl,yu,zl,zu] = deal( yl + alpha_bar * delta_yl, yu + alpha_bar * delta_yu, zl + alpha_bar * delta_zl, zu + alpha_bar * delta_zu );
        [sl,su,tl,tu] = deal( sl + alpha_bar * delta_sl, su + alpha_bar * delta_su, tl + alpha_bar * delta_tl, tu + alpha_bar * delta_tu );

    end
    disp(max(yl))
end