function [U_p, U_n]= refpotantial (theta_p, theta_n)


%% Data fit

 U_p = 3.4323 - 0.8428.*exp(-80.2493.*((1-theta_p).^1.3198)) - ...
             (3.2474*10^-6).*exp(20.2645.*((1-theta_p).^3.8003))+ ...
             (3.2482*10^-6).*exp(20.2646.*((1-theta_p).^3.7995));   
 
 U_n=0.6379 + 0.5416.*exp(-305.5309*theta_n) + ...
            0.044.*tanh(-(theta_n-0.1958)/0.1088) - ...
            0.1978.*tanh((theta_n-1.0571)/0.0854) - ...
            0.6875.*tanh((theta_n+0.0117)/0.0529) -...
            0.0175.*tanh((theta_n-0.5692)/0.0875);

    
       
end