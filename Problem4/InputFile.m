close all; clear all; 
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

function [c,dc,ddc] = EQcon(x);
    tmp = (x(1)+2);
    c = tmp^2-x(2)-3;
    if nargout ==2
        dc = [2*tmp;-1];
        ddc = [];
    elseif nargout ==3
        dc = [2*tmp;-1];
        ddc = zeros(2,2,1);
        ddc(1,1,1) = 2;
    else
        dc = [];
        ddc = [];
    end
end


function [c,dc,ddc] = InEQcon(x);
    c = [-4*x(1)+10*x(2)];
    if nargout ==2
        dc = [-4; 10];
        ddc = [];
    elseif nargout ==3
        dc = [-4; 10];
        ddc = zeros(2,2,1);
    else
        dc = [];
        ddc = [];
    end
end

% Initial guess:
x0 = [-3; 0]; y0 = 1;
xl = [-4;-2.5]; xu = [1;4];
cl = [-18]; cu = [34];
contourplt(@objective,@EQcon,@InEQcon,cl,cu,xl,xu,[-4.5,4.5,-4.5,4.5],100);
exportgraphics(gcf,'Himmelblau.png','Resolution',300)
[x] = fminconWrapper(x0);
exportgraphics(gcf,'Himmelblau_fmincon_-3_0.png','Resolution',300)
contourplt(@objective,@EQcon,@InEQcon,cl,cu,xl,xu,[-4.5,4.5,-4.5,4.5],100);
[x] = IpoptWrapper(x0);
exportgraphics(gcf,'Himmelblau_ipopt_-3_0.png','Resolution',300)
contourplt(@objective,@EQcon,@InEQcon,cl,cu,xl,xu,[-4.5,4.5,-4.5,4.5],100);
[x] = SQPLineSearchSolver(x0,@objective,@EQcon,@InEQcon,cl,cu,xl,xu,false);
exportgraphics(gcf,'Himmelblau_linesearch_-3_0.png','Resolution',300)
contourplt(@objective,@EQcon,@InEQcon,cl,cu,xl,xu,[-4.5,4.5,-4.5,4.5],100);
[x] = SQPLineSearchSolver(x0,@objective,@EQcon,@InEQcon,cl,cu,xl,xu,true);
exportgraphics(gcf,'Himmelblau_linesearch_BFGS_-3_0.png','Resolution',300)

