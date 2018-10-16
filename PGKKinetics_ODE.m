function [NADH] = PGKKinetics_ODE( k,x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   k = [490 860 1100 53]; rates from Makhatadze plos paper
% ATP + 3PG + PGK <==> PGK + ADP + 3BPG; rate = k(1), k(2) (kforward,
% kreverse)
% 3BPG + NADH + GAPDH <==> GAP + NAD+ + GAPDH; rate = k(3), k(4) (kforward,
% kreverse)

    A0 = 5E-3; %concentration of 3PG
    B0 = 5E-3; %concentration of ATP
    C0 = 0.1E-3; %concentration of NADH
    PGK = 5E-9;
    GAPDH = 1E-6;
    
    [t,y] = ode45(@PGKKinetics_eqns, x, [0 0]);
    NADH = C0 - y(:,2);
    
%     figure
%     plot(t, A0 - y(:,1),'b');
%     hold on
%     plot(t, B0 - y(:,1),'r');
%     legend('3PG', 'ATP');
%     legend('boxoff');
%     
%     figure
%     plot(t, C0 - y(:,2),'k');
%     legend('NADH');
%     legend('boxoff');

        function dydt = PGKKinetics_eqns( t, y )
        %UNTITLED Summary of this function goes here
        %   Detailed explanation goes here
        %

            dydt = [(k(1)*(A0-y(1))*(B0-y(1))*PGK) - (k(2)*(y(1)-y(2))*y(1)*PGK); (k(3)*(y(1)-y(2))*(C0-y(2))*GAPDH) - (k(4)*GAPDH*y(2)^2)];
        end

end

