function KGK_msd(D,Deff,L,t,x0,N)

%-------------------------------------------------------------------------------------
% This file contains the Matlab code mentioned in:
% "Molecular Motion in Cell Membranes: Analytic Study of Fence-Hindered Random Walks"
% Vasudev M. Kenkre, Luca Giuggioli, and Ziya Kalay, Phys. Rev. E 77, 051907 (2008)
% 
% 
% Parameters are:
% 
% D    : short time diffusion constant (as in the paper)
% Deff : effective diffusion constant (as in the paper)
% L    : compartment size (as in the paper)
% t    : time vector
% x0   : initial position of the random walker ( the first compartment is defined by
%        -L/2 < x < L/2 )
% N    : 2N+1 is the number of poles for which the residues are going to be calculated
%        in order to perform the inverse Laplace transform (precision of finding roots
%        should be modified if N > 1000, see the function 'findroots' below)
%
% Example:
%
% KGK_msd(3,0.2,0.1,logspace(-6,2,100),0.04,1000);
%-------------------------------------------------------------------------------------


d = Deff/(D-Deff); % d is the lowercase gamma in the paper (d = 0 --> complete confinement)
a = 2*x0/L; % a is lowercase alpha in the paper


% When comparing theory to experiments, uncomment #1 below (and comment %2 )and
% take 't' as a vector whose elements are real time points as in the
% experiment.
% 
% And when using dimensionless time for all other purposes uncomment #2 (and comment #1). 

% #1 : when comparing the theory to experiments
tt=4*D*t/L/L;
q=1/2;

% #2 : when using dimensionless time :
% tt=t;
% q=2/L/L;

%-------------------

if d ~= 0;
    rts=findroots(1/d,N); % N is the number of nonzero roots. Roots of the tanh equation 
    % (see the function findroots(x,y) at the end of this file)
    
    s2=i*rts(2:(N+1));
    s1=i*pi*(1:N); % roots of sinh(s)
elseif d == 0;
    s2=i*pi*0.5*(2*(0:N-1)+1); % roots of cosh(s)
    s1=i*pi*(1:N); % roots of sinh(s)
end;


% average MSD
% tt is the dimensionless time, 4*D/L^2*t

if d~=0
    res0=0.5*(2*d+6*tt*d+6*tt)/(3+3*d^2+6*d); % residue at s=0 for all times
elseif d==0
    res0=(1/6)*(-2+6*tt);
end

res=zeros(1,length(tt));

if d~=0
    n=1;
    for T=tt;
        res(n)=sum(exp(s2.*s2*T)./(s2.*s2.*(1+d*(1-tanh(s2).^2))));
        n=n+1;
    end;
    res=2*real(res);
    MSDavg=L*L*0.5*(tt - 1/d*(tt - res0 - res ));

elseif d==0;
    n=1;
    for T=tt;
        res(n)=sum(exp(s2.*s2*T)./(s2.*s2.*s2.*s2));
        n=n+1;
    end;
    res=2*real(res);
    MSDavg=L*L*0.5*(tt - (res0 + res));

end;

%-----------------------------------------------------------------

% MSD when starting from a particular initial position x0, denoted by MSDx0
% throughout. 

res=zeros(1,length(tt));

res01 = a*a*0.5 - 1/6;
res02 = d*((6*tt*(d+1) + 3*a*a*(d+1) - (d+3))/(6*(d+1)^2));
res03 = a*a/(d+1);

% Below, T is a scalar related to dimensionless time

if d ~= 0;

    n=1;
    for T=tt;
        res(n) = sum([-(cosh(a*s1).*exp(s1.*s1*T)./(s1.*s1.*cosh(s1))) ...
            ( exp(s2.*s2*T).*( s2*a.*sinh(a*s2) + d*cosh(a*s2) ) ...
            ./ (s2.*s2.*( s2.*sinh(s2) + (1+d)*cosh(s2) ) ) )]);
        n=n+1;
    end;

elseif d==0;

    n=1;
    for T=tt;

        res(n) = sum(-(cosh(a*s1).*exp(s1.*s1*T)./(s1.*s1.*cosh(s1))) ...
            + a*(sinh(a*s2).*exp(s2.*s2*T))./(s2.*(cosh(s2)+s2.*sinh(s2)) ));

        n=n+1;
    end;
end;

res=2*real(res); % because half of the poles were considered
MSDx0= L*L*0.5*( - ( res01 - res02 - res03 - res ) );%- res3 ) );

% for comparison purposes, the routine below evaluates the analytic solution for the case
% of complete confinement. Below, 'excc' is the MSDx0 for complete
% confinement calculated exactly (this routine can be skipped without any problems)

if d==0;

    % exact solution for complete confinement

    excc=zeros(1,length(tt));
    for m=0:1000; % 1000 terms are taken in the summation, increase for accuracy.
        excc=excc - 8*a*(-1)^m*sin(a*(2*m+1)*pi*0.5).*exp(-pi*pi*0.25*(2*m+1)^2*tt)/(pi*pi*(2*m+1)^2) ...
            + 2*(-1)^(m+1)*cos(a*(m+1)*pi)*exp(-(m+1)^2*pi*pi*tt)/((m+1)^2*pi*pi)    ;
    end;

    excc= L*L*0.5*(excc + (1 + 3*a*a)/6 );

end


%----------------------
% D(tau)/D computed according to Eq. (31) in the paper 
% This is for a specified initial condition x0. One should be able to
% average Eq. (31) over the initial conditions, contained in 'a' below, to get
% D*tau)/D average in a closed form. 
%----------------------

res = zeros(1,length(tt));

% the last term in D(tau)/D in the paper
n=1;
for T=tt;
    res(n) = sum( exp(real(s2.*s2*T)) .* ( a*s2.*sinh((a*s2)) + d*cosh((a*s2)) )...
        ./ (s2.*sinh((s2)) + (d+1)*cosh((s2)) )   );
    n=n+1;
end;

D_tau_over_D = ( d/(d+1) + 1 - t_4(a*0.5*pi,i*pi*tt,10000)  + 2*real(res) ); 
% 2 times res because half of the poles were considered and the other half is the same.
% Here we took 10000 terms in the summation for the theta_4 function,
% increase for accuracy.


%----------------------
% D(tau)/D (Average) computed according to Eq. (31) in the paper 
% This is for average over all initial conditions. 
%----------------------

res = zeros(1,length(tt));

n=1;
for T=tt;
    res(n) = sum( exp(real(s2.*s2*T)) .* ( (cosh(s2) - sinh(s2)./s2 ) + ...
        d*sinh(s2)./s2 )  ./ (s2.*sinh((s2)) + (d+1)*cosh((s2)) )   );
    n=n+1;
end;

D_tau_over_D_avg = ( d/(d+1) + 1 - t_4_avg(i*pi*tt,10000)  + 2*real(res) ); 
% 2 times res because half of the poles were considered and the other half is the same.
% Here we took 10000 terms in the summation for the integral of theta_4 function over a,
% in [-1,1]. 

% OUTPUT --------------------

    close all;
    figure;
        
    loglog(t,q*real(MSDavg),'-k');
    xlabel('time');
    ylabel('MSD averaged over all positions');
    
    figure;
    loglog(t,q*real(MSDx0),'-k');
    xlabel('time');
    ylabel(['MSD with initial position, x0 =', num2str(x0) ]);
    
    figure;
    loglog(t,real(D_tau_over_D),'-r');
    xlabel('time');
    ylabel('D(\tau)/D');
    
    figure;
    loglog(t,real(D_tau_over_D_avg),'-r');
    xlabel('time');
    ylabel('D(\tau)_{avg}/D');
    
% END OF OUTPUT --------------------


% END OF MAIN PROGRAM --------------

%------------

function theta_4 = t_4(u,tau,N)

theta_4 = zeros(1,length(tau));

for n=-N:N;
    theta_4 = theta_4 + (-1)^n .* exp(i*pi*tau*n*n) .* exp(2*i*n*u) ;
end

%---------------
function theta_4 = t_4_avg(tau,N)

theta_4 = zeros(1,length(tau));

for n=[-N:-1:-1 1:N];
    theta_4 = theta_4 + (-1)^n .* exp(i*pi*tau*n*n) .* sin(pi*n)/(pi*n) ;
end
theta_4 = theta_4 + 1;
%---------------

function rt=findroots(a,N)

% this function finds the roots of -ay=tan(y) by bisection method
% the locations of the roots are symmetric around the origin.

d=1e-12; % a small positive number employed in order to avoid the singularity at x=pi*0.5;
prec=1e-12; % this is the precision with which the roots are found

base=pi*0.5*(1:2:2*N-1); % N+1 is the number of roots to be found

rt = zeros(1,N+1);

for j=1:length(base);

    x(1)=base(j)+d;
    x(2)=x(1)+d;

    while tan(x(2))+a*(x(2)) < 0;
        x(2)=x(2)+0.1;
    end;

    % now we know that the root is in between x(1) and x(2)
    eps=x(2)-x(1); % this is the length of the interval in which the root lies

    while eps >= prec;

        x(3)=0.5*(x(1)+x(2));
        if tan(x(3))+a*(x(3)) < 0
            x(1)=x(3);
            eps=x(2)-x(1);
        elseif tan(x(3))+a*(x(3)) > 0
            x(2)=x(3);
            eps=x(2)-x(1);
        else
            rt(j)=x(3);
            eps=0;
        end

        rt(j)=x(3);
    end;
end;

rt=[0 rt];