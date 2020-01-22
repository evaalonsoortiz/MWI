function [A_fit, T2s_fit, resnorm] = three_comp_fit(signal, echo_times, B0)


%%-------------------------------------------------------------------------
%% Set default settings according to relaxation time
%%-------------------------------------------------------------------------
if B0 == 3
    
    %-------------------------------------------------------------------------
    % 3 T values
    %-------------------------------------------------------------------------
    T2s_MW_guess = double(10e-3);
    T2s_EW_guess = double(40e-3);
    T2s_AW_guess = double(60e-3);

    ub_MW = double(25e-3);  % Nam et al. NeuroImage 2015
    ub_EW = double(150e-3);  % Nam et al. NeuroImage 2015
    ub_AW = double(150e-3);  % Nam et al. NeuroImage 2015

    lb_MW = double(3e-3);  % Nam et al. NeuroImage 2015
    lb_EW = double(25e-3); % Nam et al. NeuroImage 2015
    lb_AW = double(25e-3);  % Nam et al. NeuroImage 2015

elseif B0 == 7 
    
    %-------------------------------------------------------------------------
    % 7 T values
    %-------------------------------------------------------------------------
    T2s_MW_guess = double(6e-3);
    T2s_EW_guess = double(30e-3);
    T2s_AW_guess = double(40e-3);

    ub_MW = double(15e-3);  
    ub_EW = double(150e-3); 
    ub_AW = double(150e-3);  

    lb_MW = double(3e-3);  
    lb_EW = double(15e-3); 
    lb_AW = double(15e-3);  

    %-------------------------------------------------------------------------
end

A_MW_guess = double(0.1*abs(signal(1)));
A_EW_guess = double(0.5*abs(signal(1)));
A_AW_guess = double(0.5*abs(signal(1)));

guess = [A_MW_guess, T2s_MW_guess, A_EW_guess, T2s_EW_guess, A_AW_guess, T2s_AW_guess];

lb_A_MW = 0;
lb_A_EW = 0;
lb_A_AW = 0;

ub_A_MW = double(0.3*signal(1));
ub_A_EW = double(1.5*signal(1));
ub_A_AW = double(1.5*signal(1));

lb = [lb_A_MW, lb_MW, lb_A_EW, lb_EW, lb_A_AW, lb_AW];
ub = [ub_A_MW, ub_MW, ub_A_EW, ub_EW, ub_A_AW, ub_AW];


%%-------------------------------------------------------------------------
%% data fitting and analysis
%%-------------------------------------------------------------------------

[result, resnorm, res, flag] = lsqnonlin(@calc_diff,guess,double(lb),double(ub));

A_fit.MW = result(1);
T2s_fit.MW = result(2); 

A_fit.EW = result(3);
T2s_fit.EW = result(4); 

A_fit.AW = result(5);
T2s_fit.AW = result(6); 
            
%%-------------------------------------------------------------------------
%% function definitions
%%-------------------------------------------------------------------------

function diff = calc_diff(guess)
    
    diff = squeeze(double(signal)) - (guess(1)*exp(-echo_times/guess(2))+guess(3)*exp(-echo_times/guess(4))+guess(5)*exp(-echo_times/guess(6)));

end


end
