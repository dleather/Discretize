%% MSVAR Tauchen (1986): Uniform grid

%% Discretization Paramters
nS = 2; %Ns
N = 2; %Number of elements in VAR
nPoints = NaN(N,nS);
nPoints = [5,5;5,5]; %Number of points of variable j in regime i
nSig = 3.25; %Number of std. dev. from mean to make grid;

%% MS-VAR(1) Specification
Pi = [0.97,0.03;0.95,0.05];
PhiCell = cell(nS);
PhiCell{1} = [0.7,0.1;0,0.9];
PhiCell{2} = [0.1,0;0,0.99];

%Check for MSS - following Cho (2016, RED)
    Rss = [Pi(1,1).*kron(PhiCell{1},PhiCell{1}),...
        Pi(2,1).*kron(PhiCell{1},PhiCell{1});...
        Pi(1,2).*kron(PhiCell{2},PhiCell{2}),...
        Pi(2,2).*kron(PhiCell{2},PhiCell{2})];
    
    if max(eig(Rss))<1
        disp('MS-VAR is MSS')
    else
        error('MS-VAR specification not MSS as defined in Cho (2016)')
    end
    

