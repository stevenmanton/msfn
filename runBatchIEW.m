function iewAll = runBatchIEW(varargin)

tic

%% ----- Interpret the input flags -----

% Parse the inputs:

p = inputParser;
p.addParamValue('N',-1,@(x) isnumeric(x) && x>0);
p.addParamValue('R',2,@(x) isnumeric(x) && all(x>0));
p.addParamValue('W',1,@(x) isnumeric(x) && all(x>0));
p.addParamValue('lam',-1,@(x) isnumeric(x) && all(x>0));
p.addParamValue('b',-1,@(x) isnumeric(x) && all(x>0));
p.parse(varargin{:});
N = p.Results.N;
R = p.Results.R;
W = p.Results.W;
lam = p.Results.lam;
b = p.Results.b;


lamSpecified = all(lam ~= -1);
bSpecified = all(b ~= -1);

if N == -1
   N = max([length(R),length(W)]);
end

% Set up the geometry parameters:
% if N == 1
%    Rs = R(1);
%    Ws = W(1);
% else
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
   
   if bSpecified % only do if lam is specified
      if length(b) == 1
         bs = b*ones(1,N);
      elseif length(b) == 2
         bs = logspace(log10(b(1)),log10(b(2)),N);
      elseif length(b) > 2
         bs = b;
      end
   end
% end

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
   
   iewAll(i).cren.N = [18 18];
   % iewAll(i).cren.N = [10 10];
   iewAll(i).plotOptimizeCrenWidths = false;
%    iewAll(i).run('maxIter',3,'plot')
   iewAll(i).run('maxIter',3,'calcbsurf','iter','calcbedge','end');
%    iewAll(i).run('maxIter',3);
   
   fprintf(1,'<Phi^2> for edges: \t%g\t%g\t%g\n', iewAll(i).Phi2edge);
end

toc

end