function [x] = fminconWrapper(x0);
    options = optimoptions('fmincon', 'OutputFcn', @outfun, 'Display', 'iter');
    xl = [-4;-2.5]; xu = [1;4];
    [x] = fmincon(@objective,x0,[],[],[],[],xl,xu,@nonlcon,options);
    global iter_data;
    figure(1);
    hold on;
    plot(iter_data(:,2),iter_data(:,3),':ok')
end


function stop = outfun(x, optimValues, state)
    global iter_data;
    stop = false; % Continue the optimization

    switch state
        case 'init'
            % Initialize the storage
            iter_data = [];
        case 'iter'
            % Store the current iteration data
            iter_data = [iter_data; optimValues.iteration, x', optimValues.fval];
        case 'done'
            % Finalize the storage (if needed)
    end
 end

 function [f,df,ddf] = objective(x);
        tmp1 = (x(1)^2+x(2)-11);
        tmp2 = x(1)+x(2)^2-7;
        f = tmp1^2+tmp2^2;
        if nargout ==2
            df = [4*x(1)*tmp1+2*tmp2;
                  2*tmp1+4*x(2)*tmp2];
            ddf = [];
        elseif nargout ==3
            df = [4*x(1)*tmp1+2*tmp2;
                  2*tmp1+4*x(2)*tmp2];
            ddf = [4*tmp1+8*x(1)^2+2,4*(x(1)+x(2));
                   4*(x(1)+x(2)), 4*tmp2+8*x(2)^2+2;];
        else
            df = [];
            ddf = [];
        end
 end

 function [c,ceq] = nonlcon(x)
    c = -[-4*x(1)+10*x(2)+18;-(-4*x(1)+10*x(2))+34];
    tmp = (x(1)+2);
    ceq = tmp^2-x(2)-3;
end