function iewAll = runBatchIEW(varargin)
% This function interprets and executes a batch of simulations specified by
% batchIEW.

tic

%% ----- Interpret the input flags -----

% Parse the inputs:

p = inputParser;
p.addParamValue('N',-1,@(x) isnumeric(x) && x>0);
p.addParamValue('R',2,@(x) isnumeric(x) && all(x>0));
p.addParamValue('W',1,@(x) isnumeric(x) && all(x>0));
p.addParamValue('lam',-1,@(x) isnumeric(x) && all(x>0));
p.addParamValue('b',-1,@(x) isnumeric(x) && all(x>0));
p.addParamValue('crenN',-1,@(x) (isnumeric(x) && all(x>0)) || iscell(x));
p.addParamValue('iewArgs',{},@iscell)
p.parse(varargin{:});
N = p.Results.N;
R = p.Results.R;
W = p.Results.W;
lam = p.Results.lam;
b = p.Results.b;
crenN = p.Results.crenN;
iewArgs = p.Results.iewArgs;

% Determine if lambda and b were specified in the command line:
lamSpecified = all(lam ~= -1);
bSpecified = all(b ~= -1);

% --- Determine N from input parameters if not specified ---
if N == -1 % N not specified
   if length(R) > 1 && length(W) > 1 % Both R and W vary
      assert(length(R) == length(W),'Length(R) ~= length(W)')
      N = length(R);
   elseif length(R) > 1 % Just R varies; W is fixed
      N = length(R);
   elseif length(W) > 1 % Just W varies; R is fixed
      N = length(W);
   elseif length(crenN) > 1 % Geometry fixed, crenN varies
      assert(length(R) == 1 && length(W) == 1,'Geometry should be fixed.')
      N = length(crenN);
   else % Nothing is recognized as varying
      N = 1;
   end
end

% Set up the geometry parameters:
if length(R) == 1
   Rs = R*ones(1,N);
elseif length(R) == 2
   Rs = logspace(log10(R(1)),log10(R(2)),N);
elseif length(R) > 2
   Rs = R;
end

if length(W) == 1
   Ws = W*ones(1,N);
elseif length(W) == 2
   Ws = logspace(log10(W(1)),log10(W(2)),N);
elseif length(W) > 2
   Ws = W;
end

if lamSpecified % only do if lam is specified
   if length(lam) == 1
      lams = lam*ones(1,N);
   elseif length(lam) == 2
      lams = logspace(log10(lam(1)),log10(lam(2)),N);
   elseif length(lam) > 2
      lams = lam;
   end
end

if bSpecified % only do if b is specified
   if length(b) == 1
      bs = b*ones(1,N);
   elseif length(b) == 2
      bs = logspace(log10(b(1)),log10(b(2)),N);
   elseif length(b) > 2
      bs = b;
   end
end

% Interpret the number of crenellations crenN:
if length(crenN) == 1 && crenN == -1
   % Default crenN = [18 18]:
   crenNmat =  repmat([18 18],N,1);
elseif length(crenN) == 1 && crenN ~= -1
   % Single crenN specified, e.g. crenN = 12:
   crenNmat =  repmat(crenN*[1 1],N,1);
elseif length(crenN) >= 2
   % All crenN specified as a vector:
   crenNmat = repmat(reshape(crenN,N,1),1,2);
end

% Default parameters for iewasher:
if isempty(iewArgs)
   iewArgs = {'maxIter',3,'calcbsurf','iter','calcbedge','end'};
end

%% ----- Run the batch -----

iewAll = iewasher.empty(0,0);
disp('Starting simulations...')
for i = 1:length(Ws)
   fprintf(1,'Current iteration: %d\n',i);
   iewAll(i).geom.fixedParam = 'R';
   iewAll(i).R = Rs(i);
   iewAll(i).W = Ws(i);
   
   % Write lambda if specified:
   if lamSpecified, iewAll(i).lambda = lams(i); end
   % Write b if specified:
   if bSpecified, iewAll(i).thickness = bs(i); end
   
   iewAll(i).cren.N = crenNmat(i,:);
   iewAll(i).plotOptimizeCrenWidths = false;
   %    iewAll(i).run('maxIter',3,'plot')
   % iewAll(i).run('maxIter',3,'calcbsurf','iter','calcbedge','end');
   iewAll(i).run(iewArgs{:});
   
   fprintf(1,'<Phi^2> for edges: \t%g\t%g\t%g\n', iewAll(i).Phi2edge);
end

toc

end