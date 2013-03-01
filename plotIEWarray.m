function plotIEWarray(iewAll,varargin)

p = inputParser;
p.CaseSensitive = false;
p.addParamValue('MarkerSize',6,@isnumeric);
p.addParamValue('R','out', @(x) any(strcmpi(x,{'ins','out','mid'})));
p.addParamValue('edgeComps',false);
p.addParamValue('newFig',true);

p.parse(varargin{:});

MarkerSize = p.Results.MarkerSize;
Rstr = ['R',p.Results.R];
edgeComps = p.Results.edgeComps;
newFig = p.Results.newFig;

%% Determine type of simulation in order to plot correctly:
R = [iewAll.R]'; R = round(1e3*R)/1e3;
W = [iewAll.W]'; W = round(1e3*W)/1e3;
lam = [iewAll.lambda]';
b = [iewAll.thickness]';

if length(unique(lam)) > 1
   simType = sprintf('varyLam, W = %g {\\mu}m, R = %g {\\mu}m', ...
      mean(W), mean(R));
   xToPlot = @(iew) [iew.lambda]'*1e3;
   xLab = '\lambda (nm)';
elseif length(unique(b)) > 1
   simType = sprintf('varyb, W = %g {\\mu}m, R = %g {\\mu}m', ...
      mean(W), mean(R));
   xToPlot = @(iew) [iew.thickness]'*1e3;
   xLab = 'thickness (nm)';
elseif length(unique(W)) == 1
   simType = sprintf('constW = %g {\\mu}m', mean(W));
   xToPlot = @(iew) [iew.(Rstr)]';
   xLab = [Rstr, ' ({\mu}m)'];
elseif length(unique(R)) == 1
   simType = sprintf('constR = %g {\\mu}m', mean(R));
   xToPlot = @(iew) [iew.(Rstr)]'./[iew.W]';
   xLab = 'R/W';
   % xToPlot = @(iew) 1./[iew.W]';
   % xLab = '1/W ({\mu}m^{-1})';
elseif std(R./W)/mean(R./W) < 1e-2
   simType = sprintf('constRW = %g', mean(R./W));
   xToPlot = @(iew) [iew.W]';
   xLab = 'W ({\mu}m)';
elseif std(R-W)/mean(R-W) < 1e-2
   simType = sprintf('constH = %g {\\mu}m', 2*mean(R-W));
%    xToPlot = @(iew) [iew.(Rstr)]' - mean(R-W);
%    xLab = sprintf('%s - %g {\\mu}m ({\\mu}m)',Rstr, mean(R-W));
   xToPlot = @(iew) [iew.W]';
   xLab = 'W ({\mu}m)';
else
   disp('Simulation type not recognized.')
end

%% Extract relevant quantities:

Phi2surf = reshape([iewAll.Phi2surf], 3, length(iewAll))';
Phi2edge = reshape([iewAll.Phi2edge], 3, length(iewAll))';
Phi2insEdge = reshape([iewAll.Phi2insEdge], 3, length(iewAll))';
Phi2outEdge = reshape([iewAll.Phi2outEdge], 3, length(iewAll))';
Phi2tot  = reshape([iewAll.Phi2tot] , 3, length(iewAll))'; %#ok<*UDIM>

% Check for xy surface asymmetry:
if any(abs(1-Phi2surf(:,1)./Phi2surf(:,2)) > 0.01)
   warning('batchIEW:asymmetryInPhi2xy', ...
      'The x and y-components in <Phi_surf^2> are asymmetric.')
end

%% Plot relevant traces:

leg = cell(0);
if newFig, figure, end

% % All components of MSF_surface
% loglog(xToPlot(iewAll), Phi2surf, 'o'), hold on
% leg = {leg, '\langle\Phi_x^2\rangle','\langle\Phi_y^2\rangle', ...
%    '\langle\Phi_z^2\rangle'};

% Total MSF
loglog(xToPlot(iewAll), sum(Phi2tot,2), 'k^','MarkerFaceColor','k', ...
   'MarkerSize',MarkerSize)
hold on
leg = [leg, '\langle\Phi_{tot}^2\rangle'];

% Total MSF_surface
loglog(xToPlot(iewAll), sum(Phi2surf,2), 'rs','MarkerFaceColor','r', ...
   'MarkerSize',MarkerSize)
leg = [leg, '\langle\Phi_{surf}^2\rangle'];

if ~edgeComps
   % Parallel MSF_surface
   loglog(xToPlot(iewAll), sum(Phi2surf(:,1:2),2), '*', ...
      'MarkerSize',MarkerSize)
   hold on
   leg = [leg, '\langle\Phi_{||}^2\rangle'];
end

% Total MSF_edge
loglog(xToPlot(iewAll), sum(Phi2edge,2), 'dc','MarkerFaceColor','c', ...
   'MarkerSize',MarkerSize)
leg = [leg, '\langle\Phi_{edge}^2\rangle'];

if edgeComps
   loglog(xToPlot(iewAll), sum(Phi2insEdge,2), 'p', ...
      'MarkerFaceColor',[0.478 0.063 0.894], 'MarkerSize',MarkerSize, ...
      'MarkerEdgeColor','none')
   leg = [leg, '\langle\Phi_{edge, ins}^2\rangle'];
   
   loglog(xToPlot(iewAll), sum(Phi2outEdge,2), 'o', ...
      'MarkerFaceColor',[0 0.498 0], 'MarkerSize',MarkerSize, ...
      'MarkerEdgeColor','none')
   leg = [leg, '\langle\Phi_{edge, out}^2\rangle'];
end

if ~edgeComps
   % z-component of MSF_surface
   loglog(xToPlot(iewAll), Phi2surf(:,3), 'mv', 'MarkerFaceColor','m', ...
      'MarkerSize',MarkerSize)
   leg = [leg, '\langle\Phi_z^2\rangle'];
end

if ~edgeComps
   % Analytic parallel MSF_surface:
   loglog(xToPlot(iewAll), [iewAll.analyticFun('R','ins')]','--')
   leg = [leg, '\langle\Phi_{||}^2\rangle (Analytic)'];
   
   loglog(xToPlot(iewAll), [iewAll.analyticFun('R','out')]','-.')
   leg = [leg, '\langle\Phi_{||}^2\rangle (Analytic)'];
   
   loglog(xToPlot(iewAll), [iewAll.analyticFun('R','mid')]')
   leg = [leg, '\langle\Phi_{||}^2\rangle (Analytic)'];
end

legend(leg,'location','best');


title(simType)
xlabel(xLab)
ylabel('\langle\Phi^2\rangle (\Phi_0^2)')


end