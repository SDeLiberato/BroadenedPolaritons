% Companion code for the publication: 
% Exact solution of polaritonic systems with arbitrary light and matter frequency-dependent losses
% The Journal of Chemical Physics 156, 084106 (2022)
% Erika Cortese and Simone De Liberato
% School of Physics and Astronomy, University of Southampton, Southampton, SO17 1BJ, United Kingdom

% This code calculates and plots the photonic and matter components of polaritons with
% Lorentzian broadening in both the photonic and matter resonances and a coloured photonic reservoir


clear all
close all

Coloured=0; % if different from zero considers the presence of a coloured reservoir (Sec. VI of the main paper)
Parallel=0; % if different from zero runs over a parallel pool

 if Parallel
     poolobj = gcp();
 end

%% Parameters and functions definitions

wx=1; % matter excitation frequency, used as units of all the other frequencies
g=0.3*wx; % resonant light-matter couoling
wxt=sqrt(wx^2+4*g^2); % matter frequency renormalised by diamagnetic P^2 term
zp=10^-3*wx; % numerical delta used to calculate the principal part numerically
qtol=10^-3; % tollerance for quadv

% loss rates
gammaP=0.05*wx; % photonic losses
gammaM=0.05*wx; % matter losses

% vector of photonic frequencies (k-points)
Nk=100; % number of photonic frequency considered
wk_min=0*wx; % minimal photonic frequency  
wk_max=2*wx; % maximal photonic frequency
wk_v=linspace(wk_min, wk_max, Nk); % vector of photonic frequencies 

% polariton dispersion from Eq. 8
wLP=sqrt(wk_v.^2 + wxt^2 - sqrt(wk_v.^4 + 16*g^2*wk_v.^2-2*wk_v.^2*wxt^2 + wxt^4))/sqrt(2);
wUP=sqrt(wk_v.^2 + wxt^2 + sqrt(wk_v.^4 + 16*g^2*wk_v.^2-2*wk_v.^2*wxt^2 + wxt^4))/sqrt(2);

% frequency range
Nw=200; % number of frequencies    
w_min=0*wx; % minimal frequency
w_max=2*wk_max; % maximal frequency
w_v=linspace(w_min,w_max,Nw); % frequency vector
dw=w_v(2)-w_v(1); % frequency delta

% Functions for Lorentzian resonances from Eq. 52 
eta2= @(w) 2*gammaM/pi*w./((w.^2-wxt^2).^2+gammaM^2*w.^2); % |\eta|^2
zeta2=@(w,wk) 2*gammaP/pi*w.^3./((w.^2-wk^2).^2+gammaP^2*w.^2);  % |\zeta_k|^2

% Functions for Lorentzian resonances  from Eq. 53
Z=@(w) 2*((w.^2-wxt^2))./((w.^2-wxt^2).^2+gammaM^2*w.^2) ; % Z
W=@(w,wk) 2*(wk^2*(w.^2-wk^2)-gammaP^2*w.^2)./((w.^2-wk^2).^2+gammaP^2*w.^2) ; % W_k


if Coloured

% define the parameters of the coloured reservoir 
wc=2*wx; % central frequency of the cloured resevoir
Delta=0.5*wx; % width of the cloured resevoir
kappa=0.05*wx; % coupling of the cloured resevoir
wl=wc-Delta/2; % lower band edge
wu=wc+Delta/2; % upper band edge

% specific form from Eq. 62
F=@(w) 1/Delta*heaviside(Delta/2-abs(w-wc)); 

% effective central frequency and effective losses from Eq. 61
Omega2=@(w,wk) wk^2*(1+kappa*(1-1/2*real(quadv(@(wp) wp.*F(abs(wp))./(wp-w+1i*zp),wl,wu,qtol)))); % |\Omega_k|^2
Gamma=@(w,wk) gammaP+pi/2*kappa*wk^2*F(w); 

% photonic bath functions with reservoir from Eq. 60
zeta2=@(w,wk) w.^3./((w.^2-Omega2(w,wk)).^2+w.^2.*Gamma(w,wk).^2).*(2*gammaP/pi+kappa*wk^2*F(w));

% photonic bath functions (numerical integral required)
W=@(w,wk) real(quadv(@(wp) 2*wp./(w.^2-wp.^2+1i*zp).*abs(zeta2(wp,wk)),0, w_max,qtol));

end


%% Main cycle 


% initializations
K2=zeros(Nk,Nw);
J2=zeros(Nk,Nw);

parfor ck=1:Nk % cycling over the photonic frequency
    
wk=wk_v(ck); % photonic frequency

%Coupled photonic and matter field component
K2(ck,:)=(zeta2(w_v,wk)+g^2.*eta2(w_v).*(W(w_v,wk).^2+pi^2.*zeta2(w_v,wk).^2))./((1-g^2*W(w_v,wk).*Z(w_v)).^2+...
     g^4*pi^2*(W(w_v,wk).^2.*eta2(w_v).^2+zeta2(w_v,wk).^2.*Z(w_v).^2)+g^2*pi^2*eta2(w_v).*zeta2(w_v,wk).*(2+g^2*pi^2*eta2(w_v).*zeta2(w_v,wk)));
 
J2(ck,:)=(eta2(w_v)+g^2*(Z(w_v).^2+pi^2.*eta2(w_v).^2).*zeta2(w_v,wk))./((1-g^2*W(w_v,wk).*Z(w_v)).^2+...
     g^4*pi^2*(W(w_v,wk).^2.*eta2(w_v).^2+zeta2(w_v,wk).^2.*Z(w_v).^2)+g^4*pi^2*eta2(w_v).*zeta2(w_v,wk).*(2+g^2*pi^2*eta2(w_v).*zeta2(w_v,wk)));
end

if Parallel
delete(poolobj); 
end

%% Plotting

figure(1)
mx1=max(max(K2));
pcolor(wk_v,w_v,K2'/mx1) 
shading interp
colormap(1-gray)
c=colorbar;
hold on
plot(wk_v,wk_v,'r--','Linewidth',1) % plot photonic frequency
plot(wk_v,wxt*ones(1,Nk),'b--','Linewidth',1) % plot matter frequency

if Coloured
plot(wk_v,wl*ones(1,Nk),'g:','Linewidth',1.5)
plot(wk_v,wu*ones(1,Nk),'g:','Linewidth',1.5)
end 

set(gca,'FontSize',16)
title('K_{k}^2(\omega)','FontSize',16)
xlabel('Cavity Frequency (in units of \omega_x)','FontSize',16)
ylabel('Frequency (in units of \omega_x)','FontSize',16)


figure(2)
mx2=max(max(J2));
pcolor(wk_v,w_v,J2'/mx2) 
shading interp
colormap(1-gray)
c=colorbar;
hold on
plot(wk_v,wk_v,'r--','Linewidth',1) % plot photonic frequency
plot(wk_v,wxt*ones(1,Nk),'b--','Linewidth',1) % plot matter frequency

if Coloured
plot(wk_v,wl*ones(1,Nk),'g:','Linewidth',1.5)
plot(wk_v,wu*ones(1,Nk),'g:','Linewidth',1.5)
end 

set(gca,'FontSize',16)
title('J_{k}^2(\omega)','FontSize',16)
xlabel('Cavity Frequency (in units of \omega_x)','FontSize',16)
ylabel('Frequency (in units of \omega_x)','FontSize',16)


