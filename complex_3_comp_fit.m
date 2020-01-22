function [RE_fit_multi, IM_fit_multi] = complex_3_comp_fit(signal_RE, signal_IM, echo_times, B0)

%%-------------------------------------------------------------------------
%% Set default settings according to relaxation time
%%-------------------------------------------------------------------------

signal_mag = sqrt(signal_RE.^2+signal_IM.^2);
[A_fit, T2s_fit, resnorm] = three_comp_fit(signal_mag, echo_times,B0);

if B0 == 3

    T2s_MW_guess = double(T2s_fit.MW);
    T2s_EW_guess = double(T2s_fit.EW);
    T2s_AW_guess = double(T2s_fit.AW);
    
%     omega_MW_guess = 0;
%     omega_AW_guess = 0;

    omega_MW_guess = 2*pi*12;
    omega_AW_guess = -2*pi;
    
    
    lb_MW = 3e-3;  % Nam et al. NeuroImage 2015
    lb_EW = 25e-3; % Nam et al. NeuroImage 2015
    lb_AW = 25e-3;  % Hwang et al. NeuroImage 2010, however this component was considered to be EW
    
    lb_omega_MW = 0;
    lb_omega_AW = -2*pi*5;

    ub_MW = 25e-3;  % Nam et al. NeuroImage 2015
    ub_EW = 150e-3;  % Hwang et al. NeuroImage 2010, however this component was considered to be EW
    ub_AW = 150e-3;  % Nam et al. NeuroImage 2015 
    
    ub_omega_MW = 2*pi*12;
    ub_omega_AW = 0;

    ub_A_MW = 0.3*abs(signal_RE(1)+1i*signal_IM(1));
    ub_A_EW = 1.5*abs(signal_RE(1)+1i*signal_IM(1));
    ub_A_AW = 1.5*abs(signal_RE(1)+1i*signal_IM(1));
    
elseif B0 == 7 

    T2s_MW_guess = double(T2s_fit.MW);
    T2s_EW_guess = double(T2s_fit.EW);
    T2s_AW_guess = double(T2s_fit.AW);
    
    omega_MW_guess = 0;
    omega_AW_guess = 0;

%     omega_MW_guess = 2*pi*25;
%     omega_AW_guess = -2*pi;
    
    lb_MW = 3e-3;   
    lb_EW = 15e-3;  
    lb_AW = 15e-3;  
    
    lb_omega_MW = 0;
    lb_omega_AW = -2*pi*8;

    ub_MW = 15e-3;  % Nam et al. NeuroImage 2015
    ub_EW = 150e-3;  % Hwang et al. NeuroImage 2010, however this component was considered to be EW
    ub_AW = 150e-3;  % Nam et al. NeuroImage 2015
    
    ub_omega_MW = 2*pi*25;
    ub_omega_AW = 0;
    
    ub_A_MW = 0.3*abs(signal_RE(1)+1i*signal_IM(1));
    ub_A_EW = 1.5*abs(signal_RE(1)+1i*signal_IM(1));
    ub_A_AW = 1.5*abs(signal_RE(1)+1i*signal_IM(1));

end

A_MW_guess = A_fit.MW;
A_EW_guess = A_fit.EW;
A_AW_guess = A_fit.AW;

guess_RE = [A_MW_guess, T2s_MW_guess, omega_MW_guess, A_EW_guess, T2s_EW_guess, A_AW_guess, T2s_AW_guess, omega_AW_guess];
guess_IM = [A_MW_guess, T2s_MW_guess, omega_MW_guess, A_AW_guess, T2s_AW_guess, omega_AW_guess];

lb_A_MW = 0;
lb_A_EW = 0;
lb_A_AW = 0;

lb_RE = [lb_A_MW, lb_MW, lb_omega_MW, lb_A_EW, lb_EW, lb_A_AW, lb_AW, lb_omega_AW];
ub_RE = [ub_A_MW, ub_MW, ub_omega_MW, ub_A_EW, ub_EW, ub_A_AW, ub_AW, ub_omega_AW];

lb_IM = [lb_A_MW, lb_MW, lb_omega_MW, lb_A_AW, lb_AW, lb_omega_AW];
ub_IM = [ub_A_MW, ub_MW, ub_omega_MW, ub_A_AW, ub_AW, ub_omega_AW];


%%-------------------------------------------------------------------------
%% data fitting and analysis
%%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Use lsqcurvefit with MultiStart
%--------------------------------------------------------------------------

options=optimset('TolX',1e-5,'TolFun',1e-5);

fitfcn_RE = @(x,xdata)(x(1)*exp(-xdata./x(2)).*cos(xdata.*x(3))+x(4)*exp(-xdata/x(5))+x(6)*exp(-xdata./x(7)).*cos(xdata.*x(8)));
problem = createOptimProblem('lsqcurvefit','x0',double(guess_RE),'objective',fitfcn_RE,'lb',double(lb_RE),'ub',double(ub_RE),'xdata',echo_times,'ydata',squeeze(double(signal_RE)));%,'options',options);
ms = MultiStart('UseParallel',true);
[RE_fit_multi,err_RE_fit_multi] = run(ms,problem,2000);

fitfcn_IM = @(x,xdata)(-x(1)*exp(-xdata./x(2)).*sin(xdata.*x(3))-x(4)*exp(-xdata./x(5)).*sin(xdata.*x(6)));
problem = createOptimProblem('lsqcurvefit','x0',double(guess_IM),'objective',fitfcn_IM,'lb',double(lb_IM),'ub',double(ub_IM),'xdata',echo_times,'ydata',squeeze(double(signal_IM)));%,'options',options);
ms = MultiStart('UseParallel',true);
[IM_fit_multi,err_IM_fit_multi] = run(ms,problem,2000);



end
