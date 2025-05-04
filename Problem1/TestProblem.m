clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 1: Equality Constrained Convex QP (subproblem 6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(42069); % Set seed for reproducability
npoint = 40;
ftsize = 12; % Font size
bs = round(linspace(8.5,18.68,npoint));  % Problem sizes
solutionX = zeros(length(bs),5);
solutionLambda = zeros(length(bs),2);
solutionValue = zeros(length(bs),1);

H =[6.0000 1.8600 1.2400 1.4800 -0.4600
    1.8600 4.0000 0.4400 1.1200 0.5200
    1.2400 0.4400 3.8000 1.5600 -0.5400
    1.4800 1.1200 1.5600 7.2000 -1.1200
    -0.4600 0.5200 -0.5400 -1.1200 7.8000];
g = [   -16.1000
        -8.5000
        -15.7000
        -10.0200
        -18.6800];
A = [   16.1000 1.0000
        8.5000 1.0000
        15.7000 1.0000
        10.0200 1.0000
        18.6800 1.0000];
b = [15;
     1 ];

[x1,lambda1] = EqualityQPSolver(H,g,A,b,'LUdense');
disp('Optimal solution:')
disp(x1')
disp('With the optimial value')
disp(1/2*x1'*H*x1+g'*x1)

for i = 1:length(bs)
    b(1) = bs(i);
    [x1,lambda1] = EqualityQPSolver(H,g,A,b,'LUdense');
    solutionX(i,:) = x1;
    solutionLambda(i,:) = lambda1;
    solutionValue(i) = 1/2*x1'*H*x1+g'*x1;
end

figure(1)
plot(bs,solutionX)
xlabel("$b(1)$",'FontSize',ftsize,'Interpreter','latex')
ylabel("$x$ value",'FontSize',ftsize,'Interpreter','latex')
yyaxis right
plot(bs,solutionValue,'--k')
ylabel("Objective value $\phi$",'FontSize',ftsize,'Interpreter','latex')
legend(["$x(1)$","$x(2)$","$x(3)$","$x(4)$","$x(5)$","$\phi$"],'location','southoutside','NumColumns',6,'FontSize',ftsize,'Interpreter','latex')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
exportgraphics(gcf,'Problem1_6.png','Resolution',300)
