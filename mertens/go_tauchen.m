initscript % can replace with clear all, close all etc
initwrap   % please skip if not latexwrapper installed, same for finishwrap at the end

%% setup
if false
   load
else
   load aqueductdataQ.mat
%    vardata = [diff(Y(:,[1 2])) Y(2:end,7:8)]; 
   vardata = diff(Y(:,[1 2])); 
   VAR0 = VARls(vardata,1);
   F     = VAR0.F;
   F0    = VAR0.F0;
   Sigma = VAR0.Omega;
end

plotY(Y, dates, 10, Ylabels)

Ns        = 6;
Ny        = size(vardata, 2);
bandwidth = 3;
profile off

%% main
[P, y, x, MVndx, Moments0] = tauchen(F, F0, Sigma, Ns, bandwidth);

% [P, y, x, MVndx, Moments0, p0] = tauchen(F, F0, Sigma, Ns, bandwidth);

Nstar = length(P);
p0    = mclimit(P);


%% some output
dispcalc('p0'' * x')  
dispcalc('Moments0.EX''')

dispcalc('mcvar(P,x,p0)')  
dispcalc('Moments0.VarX''')

dispcalc('p0'' * y')  
dispcalc('Moments0.EY''')

dispcalc('mcvar(P,y,p0)')  
dispcalc('Moments0.VarY''')
dispcalc('disclyap(VAR0.F, VAR0.Omega)')
dispcalc('cov(vardata(2:end,:),1)')

%% p0
if false
   p0grid = mclimit(P);
   mae(p0grid,p0)
   plot([p0grid,p0])
   figure
   plot(p0grid, p0, '.')
   figure
   plot(cumsum([p0grid,p0]))
end

%% call statistics
ybar = mcmean(p0', y)
ybarc = mcmean(P, y);
dispcalc('p0'' * ybarc')

dispcalc('[ybar'' inv(eye(Ny)-F)*F0]')
dispcalc('mae(ybar, inv(eye(Ny)-F)*F0)')

%% variances
tic
Omega  = mcvar(P, y, p0);
Omega0 = disclyap(F, Sigma);
dispcalc('mae(Omega, Omega0)')

Omegas = repmat(NaN, [size(Omega), Nstar]);
for s = 1 : Nstar
   Omegas(:,:,s) = mcvar(P, y, s);
end
if Nstar < 10
   dispcalc('Omegas')
end
EOmegas     = sum(Omegas .* repmat(reshape(p0, [1 1 Nstar]), [Ny Ny 1]), 3);
OmegaYbars  = mcvar(P, P * y - repmat(p0' * y, Nstar, 1), p0);
dispcalc('EOmegas + OmegaYbars - Omega')
toc

%% Retrieve VAR parameter
tic
ACF0 = mcvar(P, y, p0, 1);
dispcalc('ACF0'' * inv(Omega)')
dispcalc('F')
dispcalc('mae(ACF0'' * inv(Omega), F)')
dispcalc('mae(mean(Omegas,3), Sigma)')
dispcalc('mae(median(Omegas,3), Sigma)')
dispcalc('std(Omegas,1,3)')
toc

%% simulate something
tic
states = mcdrawstates(P, 50100, p0);
Ysim = mcsim(states, y);
newfig
plotY(Ysim(100:500,:))
wrapcf('ysim', wrap)
simVAR = VARls(Ysim(100:end,:), 1);
dispcalc('simVAR.F')
dispcalc('F')
dispcalc('simVAR.Omega')
dispcalc('Sigma')
% VARprint(simVAR)
toc

%% wrap
finishwrap
finishscript
