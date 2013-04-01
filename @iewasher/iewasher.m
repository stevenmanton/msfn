classdef iewasher < handle
   %iewasher Class calculates mean square flux noise of square loop
   % Detailed explanation...
   %
   % iewasher Properties:
   % ----- Solution properties -----
   %    L - Computed inductance (pH)
   
   %{
   ----- A note on units -----

   Converting between current density J, magnetic field B, and flux Phi can
   get a bit tricky, especially because certain quantities are normalized.
   Below, I hope to clarify my conventions for units so that we can
   properly keep track of what's going on.

   Current density:

   Because the current density is normalized, we can simply say that the
   units of the current are amps/um^2, such that the total circulating
   current (integrated over the width and thickness) is 1 Amp. Obviously,
   this is not physical, but it's fine in terms of keeping track of units.

   Magnetic field: Is in units of Tesla.

   "Total" magnetic field (<B^2>*A). Is in units of Tesla*um^2

   Mean square flux: Units of flux quantum.
   %}
   
   properties
      % ----- Geometry properties -----
      
      % geom - Struct with inner and outer dimensions of loop
      geom = struct('oDim',[4 4],'iDim',[2 2],'fixedParam','W')
      
      % ----- Material properties -----
      
      lambda = 0.09        % Penetration depth (um)
      thickness = 0.15     % Film thickness of loop (um)
      spinDens = 2.5e17;   % Surface spin density (1/m^2)
      
      % ----- Solution properties -----
      
      runStartTime % Time stamp set when 'run' command is executed
      
      L        % Extracted inductance (pH)
      
      % Extracted from *.inp file:
      nodes    % struct containing information on all nodes
      edges    % struct containing information on all segments
      
      Bsurf    % Vector with magnetic field info: [x, y, z, Bx, By, Bz]
      B2surf = nan(1,3)    % Components of <B_surf^2> (T^2*um^2)
      B2insEdge = nan(1,3) % Components of <B_insedge^2> (T^2*um^2)
      B2outEdge = nan(1,3) % Components of <B_outedge^2> (T^2*um^2)
      
      % [sum(Bx^2), sum(By^2), sum(Bz^2)]
      % Total magnetic field, calculated for the lower left corner
      % and multiplied by 4.
      
      
      relJ2AnalyticNumeric % Ratio of int(J_analytic^2) / int(J_numeric^2)
      % Optimum value of Van Duzer fudge factor to match numeric J(x)
      aOptimum
      relB2AnalyticNumeric % Ratio of int(B_analytic^2) / int(B_numeric^2)
      
      % Error in the expansion of epsilon (see Mathematica workbook), that
      % is, the error introduced by not integrating entirely to the edge.
      errorFromEpsExpansion
      
      % Estimated Phi^2 using analytic current density equation (Eq. (5)
      % from Bialczak)
      Phi2estimateFromAnalyticJ
      % Estimated Phi^2 using analytic magnetic field
      Phi2estimateFromAnalyticB
      
      err  % error generated during the 'run' command
      
      % ----- Numerical properties -----
      
      runIter = 1 % Current iteration of the run
      gap = struct('rel',1,'mode','rel')
      cren = struct('H',20,'line',20*ones(7,1),'hole',20*ones(6,1), ...
         'mode','default', 'N',[18 18])
      
      iterInfo % Struct with information saved from each iteration
      
      runParallel = true % run code (magnetic field) in parallel?
      
      % ----- Folder properties -----
      
      % Folder where FastHenry and InductEx executables are stored
      foldRoot = 'C:\Dropbox\Berkeley\MATLAB\github\msfn\bin\'
%       foldRoot = 'C:\test\bin\'
      
      % ----- Plotting properties -----
      
      % Plot the optimized crenellation widths versus iteration?
      plotOptimizeCrenWidths = false
   end
   
   properties (Dependent)
      % ----- Geometry properties -----
      
      R     % R measured to outside of trace (um)
      Rins  % R measured to inside of trace (um)
      Rout  % R measured to outside of trace (um)
      Rmid  % R measured to middle of the trace (um)
      
      W     % linewidth (um)
      RW    % R/W
      oDim  % Outer dimension of loop (um)
      iDim  % Outer dimension of loop (um)
      wGap  % Gap of slit (um)
      
      
      % ----- Solution properties -----
      
      % Ratio of <Phi_data^2> / <Phi_Analytic^2>
      relPhi2AnalyticNumeric
      % Ratio of <Phi_data^2> / <Phi_Analytic^2> using the analytic result
      % adjusted by the difference difference in <B^2> across the
      % linewidth. In other words, departures from a unity value are due to
      % other things, such as corners, ambiguity in the definition of R,
      % etc.
      relPhi2AdjustedAnalyticToNumeric      
      % Summary of comparison of numeric to analytic answer
      relNumericAnalyticSummary
      
      % ----- Numerical properties -----
      
      gapMax      % Maximum allowed "gap" that InductEx will allow
      
      B2toPhi2    % Conversion factor from B^2 to Phi^2:
      
      B2edge      % <B_edge^2> * A_loop (T^2*um^2)
      Phi2surf    % <Phi_surf^2> (Phi_0^2/Hz)
      Phi2insEdge % <Phi_{inside edge}^2> (Phi_0^2/Hz)
      Phi2outEdge % <Phi_{outside edge}^2> (Phi_0^2/Hz)
      Phi2edge    % <Phi_edge^2> (Phi_0^2/Hz)
      Phi2tot     % <Phi_total^2> (Phi_0^2/Hz)
      J           % current density: [x, y, Jx, Jy] Units: um, (A/um^2)
      
      objData     % Struct containing data for nondependent properties of object
   end
   
   methods
      %% Class constructor:
      function iew = iewasher()
         % Class constructor... nothing special here.
      end % function
      
      %% copy
      function iewB = copy(iewA,iewB)
         % Perform a deep copy of iewasher object
         
         % Get metadata from iewasher class (we need a list of all the
         % variables and their properties in order to perform the copy).
         % Note that this is only a shallow copy. This routine won't do a
         % recursive copy if any of the properties in the objects are
         % handles themselves.
         % Also note that if it's important to maintain the handle and just
         % copy the data (if, for instance, an existing array of objects
         % has a handle that points to iewB, and you just want to transfer
         % the data from A to B, you can specify iewB as the second
         % argument. If iewB is not specified, then a new object will be
         % generated.
         mc = ?iewasher;
         
         % Initialize the copy if necessary:
         if nargin == 1
            iewB = iewasher;
         end
         
         % Loop over properties and copy only if the ones that aren't
         % constants or dependents
%          for i = 1:length(mc.Properties)
%             if ~mc.Properties{i}.Dependent && ...
%                   ~mc.Properties{i}.Constant
%                iewB.(mc.Properties{i}.Name) = ...
%                   iewA.(mc.Properties{i}.Name);
%             end
%          end
         
         for i = 1:length(mc.PropertyList)
            if ~mc.PropertyList(i).Dependent && ...
                  ~mc.PropertyList(i).Constant
               iewB.(mc.PropertyList(i).Name) = ...
                  iewA.(mc.PropertyList(i).Name);
            end
         end
         
      end
      
      %% writeObjData
      function writeObjData(iew, objData)
         % Write object data to various fields in iewasher object
         % (Used when performing a "deep" copy of an object)
         
         fnames = fieldnames(objData);
         for i = 1:length(fnames)
            iew.(fnames{i}) = objData.(fnames{i});
         end
      end
      
      %% Get functions (geometry)
      function val = get.R(iew)
         val = iew.geom.oDim/2;
         if val(1) == val(2), val = val(1); end
      end
      
      % Aliased functions for different definitions of R:
      function val = get.Rins(iew), val = iew.R - iew.W; end
      function val = get.Rout(iew), val = iew.R; end
      function val = get.Rmid(iew), val = iew.R - iew.W/2; end
      
      function val = get.W(iew)
         val = (iew.geom.oDim - iew.geom.iDim)/2;
         if val(1) == val(2), val = val(1); end
      end
      
      function val = get.RW(iew)
         val = iew.R ./ iew.W;
      end
      
      function val = get.oDim(iew)
         val = iew.geom.oDim;
      end
      
      function val = get.iDim(iew)
         val = iew.geom.iDim;
      end
      
      function val = get.gapMax(iew)
         % Set the GapMax parameter for InductEx. If GapMax is larger then
         % iew.wGap, then sometimes InductEx bridges the gap and the
         % simulation doesn't run correctly.
         val = min([round(1000*min(iew.oDim)/20)/1000, iew.wGap]);
      end
      
      function val = get.wGap(iew)
         % Gap in the SQUID loop (proportional to the inner dimension)
         val = (1/4) * iew.iDim(2);
      end
      
      %% Get functions (solution)
      
      % Conversion factor from B^2 to Phi^2:
      function B2toPhi2 = get.B2toPhi2(iew)
         % sigma * mu_B^2 / 3 * I^2 * Phi_0^2
         % (Note that I is normalized to equal unity)
         % Wolfram code:
         % (T^2*um^2)*(5E17/m^2)*(bohr magneton)^2 / (3*A^2*flux quantum^2)
         
         B2toPhi2 = iew.spinDens * ...
            (9.27401e-24)^2 * 1e-12 / (3*(2.067834e-15)^2);
      end            
      
      % Total mean square magnetic field along inner and outer edges:
      function B2edge = get.B2edge(iew)
         B2edge = iew.B2insEdge + iew.B2outEdge;
      end
      
      % Total mean square flux due to surface spins:
      function Phi2surf = get.Phi2surf(iew)
         if isfield(iew,'runStartTime')
            Phi2surf = iew.B2toPhi2 * iew.B2surf;
         else
            % Need to add factor of 2 because object was run before spin
            % density glitch was fixed:
            Phi2surf = 2 * iew.B2toPhi2 * iew.B2surf;
         end
      end
      
      % Total mean square flux due to inside edge spins:
      function Phi2insEdge = get.Phi2insEdge(iew)
         Phi2insEdge = iew.B2toPhi2 * iew.B2insEdge;
      end
      
      % Total mean square flux due to inside edge spins:
      function Phi2outEdge = get.Phi2outEdge(iew)
         Phi2outEdge = iew.B2toPhi2 * iew.B2outEdge;
      end
      
      % Total mean square flux due to edge spins:
      function Phi2edge = get.Phi2edge(iew)
         Phi2edge = iew.B2toPhi2 * iew.B2edge;
      end
      
      % Total mean square flux due to surface and edge spins:
      function Phi2tot = get.Phi2tot(iew)
         Phi2tot = iew.Phi2surf + iew.Phi2edge;
      end
      
      % Ratio of numerically computed <Phi^2> to the analytic <Phi^2>:
      function val = get.relPhi2AnalyticNumeric(iew)
         val = iew.analyticFun / sum(iew.Phi2surf(1:2));
      end
      
      % Ratio of the numerically computed <Phi^2> to the analytic <Phi^2>
      % that can be explained by the difference in <B^2> across the
      % linewidth. A value of unity means that the contribution of corners
      % is negligible.
      function val = get.relPhi2AdjustedAnalyticToNumeric(iew)
         val = iew.relPhi2AnalyticNumeric ./ ...
            (iew.relB2AnalyticNumeric .* (1+iew.errorFromEpsExpansion));
      end
      
      function val = get.relNumericAnalyticSummary(iew)
         val = [iew.R, iew.W, ...
            iew.relJ2AnalyticNumeric, ...
            iew.relB2AnalyticNumeric, ...
            iew.relPhi2AnalyticNumeric, ...
            iew.relPhi2AdjustedAnalyticToNumeric]';
      end
      
      % Get all non-dependent and non-constant object properties:
      function objData = get.objData(iew)
         % Get meta data for the class iewasher:
         mc = ?iewasher;
         
         % Loop over properties and add them to struct objData:
         objData = struct;
         for i = 1:length(mc.PropertyList)
            if ~mc.PropertyList(i).Dependent && ...
                  ~mc.PropertyList(i).Constant
               objData.(mc.PropertyList(i).Name) = ...
                  iew.(mc.PropertyList(i).Name);
            end
         end
      end
      
      %% Set geometry functions
      
      % ------ Geometry functions ------
      function set.R(iew,val)
         switch upper(iew.geom.fixedParam)
            case 'R' % leave iDim fixed
               iew.oDim = 2*val;
            case 'W'
               Wfixed = iew.W;
               iew.oDim = 2*val;
               iew.iDim = 2*(val-Wfixed);
            case 'RW'
               RWfixed = iew.RW;
               iew.oDim = 2*val;
               iew.iDim = 2*val*(1 - 1/RWfixed);
            otherwise
         end
      end
      
      function set.W(iew,val)
         switch upper(iew.geom.fixedParam)
            case 'R'
               iew.iDim = 2*(iew.R - val);
            case 'W' % leave iDim fixed
               iew.oDim = iew.iDim + 2*val;
            case 'RW'
               RWfixed = iew.RW;
               iew.oDim = 2*val*RWfixed;
               iew.iDim = iew.oDim - 2*val;
            otherwise
         end
      end
      
      function set.RW(iew,val)
         switch upper(iew.geom.fixedParam)
            case 'R'
               iew.iDim = 2*iew.R*(1-1/val);
            case 'W' % leave iDim fixed
               Wfixed = iew.W;
               iew.oDim = 2*val*Wfixed;
               iew.iDim = iew.oDim - 2*Wfixed;
            case 'RW'
               RWfixed = iew.RW;
               iew.oDim = 2*val*RWfixed;
               iew.iDim = iew.oDim - 2*val;
            otherwise
         end
      end
      
      function set.oDim(iew,val)
         switch length(val)
            case 1, iew.geom.oDim = val*[1 1];
            case 2, iew.geom.oDim = val;
         end
      end
      
      function set.iDim(iew,val)
         switch length(val)
            case 1, iew.geom.iDim = val*[1 1];
            case 2, iew.geom.iDim = val;
         end
      end
      
      function set.geom(iew,geomIn)
         iew.geom.oDim = round(geomIn.oDim*1000/4)*4/1000;
         iew.geom.iDim = round(geomIn.iDim*1000/4)*4/1000;
         iew.geom.fixedParam = geomIn.fixedParam;
      end
      
      %% Set functions (other)
      function set.foldRoot(iew,val)
         if exist(val,'dir')
            if ~exist([val,'fasthenry.exe'],'file')
               disp('FastHenry executable (fasthenry.exe) is missing!')
            end
         else
            disp('Root folder does not exist. Code will crash.')
         end
         iew.foldRoot = val;
      end
      
      %% Write functions
      function writeCIR(iew)
         % Write the CIR file (defines ports)
         
         fid = fopen([iew.foldRoot,'washer.cir'],'wt');
         
         fprintf(fid,'* spice netlist for InductEx4\n');
         fprintf(fid,'L1    1   2   1\n');
         fprintf(fid,'*** Modified by: Steven Anton\n');
         fprintf(fid,'* NB: Ports must be numbered sequentually\n');
         fprintf(fid,'P1     1   2\n');
         fprintf(fid,'.end\n');
         
         fclose(fid);
      end
      
      function writeCIF(iew)
         % Write CIF file (geometry)
         
         oDimCIF = iew.oDim*1000;
         iDimCIF = iew.iDim*1000;
         wGapCIF = iew.wGap*1000; % width of gap (slit)
         lWid = (oDimCIF - iDimCIF)/2;
         
         % ----- Compute necessary points for washer polygon -----
         x = [0 1 1 -1 -1 1 1 0 0 0 0 0]*oDimCIF(1)/2 + ...
            [1 0 0 0 0 0 0 1 1 -1 -1 1]*iDimCIF(1)/2;
         y = [0 0 -1 -1 1 1 0 0 0 0 0 0]*oDimCIF(2)/2 + ...
            [0 0 0 0 0 0 0 0 1 1 -1 -1]*iDimCIF(2)/2 + ...
            [-1 -1 0 0 0 0 1 1 0 0 0 0]*wGapCIF/2;
         x(13) = x(1); y(13) = y(1);
         
         xy = zeros(2*length(x),1);
         xy(1:2:end) = x;
         xy(2:2:end) = y;
         xy = round(xy); % make sure all points are integers
         
         % ----- Compute necessary points for crenellations -----
         
%          iew.updateCren();
         
         % Crenelations to refine along the width of the linewidth:
         crenW = iew.cren.line;
         crenH = iew.cren.H;
         Bmat = []; %#ok<*AGROW>
         for i = 1:2:length(crenW)
            % left side, crenellations on top:
            Bmat(end+1,:) = [crenW(i), crenH, ...
               -iDimCIF(1)/2-crenW(i)/2 - sum(crenW(1:i-1)), ...
               oDimCIF(2)/2 + crenH/2];
            % right side, crenellations on top:
            Bmat(end+1,:) = [crenW(i), crenH, ...
               iDimCIF(1)/2 + crenW(i)/2 + sum(crenW(1:i-1)), ...
               oDimCIF(2)/2 + crenH/2];
            % bottom side, crenellations on left:
            Bmat(end+1,:) = [crenH, crenW(i), ...
               -oDimCIF(1)/2 - crenH/2, ...
               -iDimCIF(2)/2-crenW(i)/2 - sum(crenW(1:i-1))];
            % top side, crenellations on left:
            Bmat(end+1,:) = [crenH, crenW(i), ...
               -oDimCIF(1)/2 - crenH/2, ...
               iDimCIF(2)/2+crenW(i)/2 + sum(crenW(1:i-1))];
         end
         
         % Crenelations to refine the spacing perpendicular to the
         % linewidth along the boundary of the hole:
         crenW = iew.cren.hole;
         crenH = iew.cren.H;
         for i = 2:2:length(crenW)
            % left side, crenellations on top:
            Bmat(end+1,:) = [crenW(i), crenH, ...
               -iDimCIF(1)/2 + crenW(i)/2 + sum(crenW(1:i-1)), ...
               oDimCIF(2)/2 + crenH/2];
            %             % right side, crenellations on top:
            %              Bmat(end+1,:) = [crenW(i), crenH, ...
            %                iDimCIF(1)/2 - crenW(i)/2 - sum(crenW(1:i-1)), ...
            %                oDimCIF(2)/2 + crenH/2];
            %             % bottom side, crenellations on left:
            %              Bmat(end+1,:) = [crenH, crenW(i), ...
            %                -oDimCIF(1)/2 - crenH/2, ...
            %                -iDimCIF(2)/2 + crenW(i)/2 + sum(crenW(1:i-1))];
            % top side, crenellations on left:
            Bmat(end+1,:) = [crenH, crenW(i), ...
               -oDimCIF(1)/2 - crenH/2, ...
               iDimCIF(2)/2 - crenW(i)/2 - sum(crenW(1:i-1))];
         end
         
         Bmat = round(Bmat); % make sure all points are integers
         
         % --------------- Write the file ---------------
         
         fid = fopen([iew.foldRoot,'washer.cif'],'wt');
         
         fprintf(fid,'(! CIF file made by santon)\n');
         fprintf(fid,'\n');
         fprintf(fid,'(LAYER 1 CIF="M1");\n');
         fprintf(fid,'(LAYER 19 CIF="TERM");\n');
         fprintf(fid,'(LAYER 30 CIF="M0");\n');
         fprintf(fid,'\n');
         fprintf(fid,'(#1 WASHERDEMO1);\n');
         fprintf(fid,'\n');
         
         fprintf(fid,'DS#1 2/2;\n');
         fprintf(fid,'9 WASHERDEMO1;\n');
         fprintf(fid,'LM1;\n');
         
         % Write the crenellations:
         for i = 1:size(Bmat,1)
            fprintf(fid,'B %d %d %d %d;\n', Bmat(i,:));
         end
         
         % Write the polygon coordinates:
         fprintf(fid,'P');
         fprintf(fid,' %d',xy);
         fprintf(fid,';\n');
         
         fprintf(fid,'LTERM;\n');
%          fprintf(fid,'B %d 0 %d %d;\n', ...
%             round([lWid(1)/4, (oDimCIF(1)-lWid(1))/2, -wGapCIF/2]));
%          fprintf(fid,'B %d 0 %d %d;\n', ...
%             round([lWid(1)/4, (oDimCIF(1)-lWid(1))/2, wGapCIF/2]));
         fprintf(fid,'B %d 0 %d %d;\n', ...
            round([lWid(1), (oDimCIF(1)-lWid(1))/2, -wGapCIF/2]));
         fprintf(fid,'B %d 0 %d %d;\n', ...
            round([lWid(1), (oDimCIF(1)-lWid(1))/2, wGapCIF/2]));
         
         % Place some labels that are necessary for InductEx. The
         % location of the text label must be somewhere along the
         % terminal, but I don't think the exact location matters, so I
         % just put it in the middle.
         fprintf(fid,'94 P1+ M1 %d %d TEXT;\n', ...
            round([(oDimCIF(1)-lWid(1))/2, wGapCIF/2]));
         % This line doesn't matter:
         % fprintf(fid,'99 18 150 0 %d %d P1+ M1;\n', ...
         %     (oDimCIF(1)-lWid(1))/2, wGapCIF/2);
         fprintf(fid,'94 P1- M1 %d %d TEXT;\n', ...
            round([(oDimCIF(1)-lWid(1))/2, -wGapCIF/2]));
         % This line doesn't matter:
         % fprintf(fid,'99 18 150 0 %d %d P1- M1;\n', ...
         %     (oDimCIF(1)-lWid(1))/2, wGapCIF/2);
         
         fprintf(fid,'DF;\n'); % end definition
         fprintf(fid,'C#1;\n'); % call definition
         fprintf(fid,'E\n'); % end file
         
         fclose(fid);
      end
      
      function writeLDF(iew)
         % Write LDF file (material properties)
         
         fid = fopen([iew.foldRoot,'UCBerkeley.ldf'],'wt');
         
         fprintf(fid,'*** Layer Definition File for UC Berkeley SQUID washers\n');
         fprintf(fid,'*** Author: Coenrad Fourie\n');
         fprintf(fid,'*** Modified by: Steven Anton\n');
         fprintf(fid,'*** Last Modification: %s\n',date);
         fprintf(fid,'*** Info: Monolayer, assume Nb\n');
         fprintf(fid,'*\n');
         fprintf(fid,'$Parameters\n');
         fprintf(fid,'* Global parameters\n');
         fprintf(fid,'Units             =  1e-6\n');
         fprintf(fid,'*  CIFUnitsPerMicron lets InductEx know the CIF coordinate scale\n');
         fprintf(fid,'*    If your layout tool puts out CIF files with 100 units per micron (default), use 100\n');
         fprintf(fid,'*    If your layout tool uses 1000 units per micron (XIC if not stripped for export), use 1000\n');
         fprintf(fid,'*    If you have no idea what this means, keep CIFUnitsPerMicron = 100\n');
         fprintf(fid,'CIFUnitsPerMicron =  1000\n');
         fprintf(fid,'gapMax            =  %g\n', iew.gapMax);
         fprintf(fid,'AbsMin            =  0.002\n');
         fprintf(fid,'ProcessHasGroundPlane = FALSE\n');
         fprintf(fid,'LastDieLayerOrder =  1\n');
         fprintf(fid,'BlankAllLayer     =  60\n');
         fprintf(fid,'BlankXLayer       =  61\n');
         fprintf(fid,'BlankYLayer       =  62\n');
         fprintf(fid,'TermLayer         =  19\n');
         fprintf(fid,'TextLayer         =  18\n');
         fprintf(fid,'Lambda            =  %g\n',iew.lambda);
         fprintf(fid,'HFilaments        =  1\n');
         fprintf(fid,'Colour            =  1\n');
         fprintf(fid,'TerminalInRange   =  0.001\n');
         fprintf(fid,'$End\n');
         fprintf(fid,'*\n');
         fprintf(fid,'* LAYERS\n');
         fprintf(fid,'** Number is GDS layer number\n');
         fprintf(fid,'** Name is layer as applied in geometry input file\n');
         fprintf(fid,'** Bias is the mask-wafer offset of an object''s border in this layer\n');
         fprintf(fid,'** Thickness is the layer thickness in microns\n');
         fprintf(fid,'** Lamba is the layer''s penetration depth in microns\n');
         fprintf(fid,'** Sigma is the layer''s bulk conductivity (only for resistive layers, and not yet used)\n');
         fprintf(fid,'** Order is the layer''s order during wafer construction - the lowest layer starts at 0, but does not need to be Ground (as in ADP)\n');
         fprintf(fid,'** Mask is the mask polarity: {1 = layer objects define material\n');
         fprintf(fid,'**                             0 = layer objects not translated to model\n');
         fprintf(fid,'**                            -1 = layer objects define cutots }\n');
         fprintf(fid,'** Filmtype is the layer material typ: {S = superconductor, N = normal conductor, I = isolator, A = auxiliary/don''t care }\n');
         fprintf(fid,'** HFilaments is the number of filaments segments are divided into over the height (overrides global HFilaments)\n');
         fprintf(fid,'** Colour is the DXF colour (for viewing purposes)\n');
         fprintf(fid,'**\n');
         fprintf(fid,'*\n');
         
         fprintf(fid,'* M0  (fake ground layer to tie up order 0)\n');
         fprintf(fid,'$Layer\n');
         fprintf(fid,'Number     =     30\n');
         fprintf(fid,'Name       =     M0\n');
         fprintf(fid,'Bias       =     0\n');
         fprintf(fid,'Thickness  =     %g\n',iew.thickness);
         fprintf(fid,'Lambda     =     %g\n',iew.lambda);
         fprintf(fid,'Order      =     0\n');
         fprintf(fid,'Mask       =     1\n');
         fprintf(fid,'Filmtype   =     S\n');
         fprintf(fid,'HFilaments =     1\n');
         fprintf(fid,'Colour     =     10\n');
         fprintf(fid,'$End\n');
         fprintf(fid,'*\n');
         
         fprintf(fid,'* M1\n');
         fprintf(fid,'$Layer\n');
         fprintf(fid,'Number     =     1\n');
         fprintf(fid,'Name       =     M1\n');
         fprintf(fid,'Bias       =     0\n');
         fprintf(fid,'Thickness  =     %g\n',iew.thickness);
         fprintf(fid,'Lambda     =     %g\n',iew.lambda);
         fprintf(fid,'Order      =     1\n');
         fprintf(fid,'Mask       =     1\n');
         fprintf(fid,'Filmtype   =     S\n');
         fprintf(fid,'HFilaments =     1\n');
         fprintf(fid,'Colour     =     10\n');
         fprintf(fid,'$End\n');
         fprintf(fid,'*\n');
         
         fprintf(fid,'* TERM\n');
         fprintf(fid,'$Layer\n');
         fprintf(fid,'Number     =     19\n');
         fprintf(fid,'Name       =     TERM\n');
         fprintf(fid,'Bias       =     0\n');
         fprintf(fid,'Thickness  =     %g\n',iew.thickness);
         fprintf(fid,'Order      =     14\n');
         fprintf(fid,'Mask       =    -4\n');
         fprintf(fid,'$End\n');
         
         fclose(fid);
      end
      
      %% run
      function varargout = run(iew,varargin)
         % run Run the iterative computation algorithm
         
         % ----- Interpret the input flags -----
         p = inputParser;
         p.CaseSensitive = false;
         p.addOptional('plot',false);
         p.addParamValue('maxIter',3,@isnumeric)
         p.addParamValue('calcBSurf','no',...
            @(x) any(strcmpi(x,{'no','end','iter'})) )
         p.addParamValue('calcBEdge','no',...
            @(x) any(strcmpi(x,{'no','end','iter'})) )
         p.parse(varargin{:});
         
         if isempty(p.Results.plot), debugPlot = true;
         else debugPlot = false; end
         
         % ----- Run the iteration -----
         iew.runStartTime = datevec(now);
         
         iew.runIter = 0;
         hitErr = false;
         while hitErr == false && iew.runIter < p.Results.maxIter
            try % Try to iterate
               % Create a temporary object while we try this iteration. If
               % it succeeds, we'll keep it. If the iteration fails, we'll
               % throw it away and return the last successful copy.
               iewTemp = copy(iew);
               
               % Write the necessary files:
               iewTemp.writeCIR();
               iewTemp.writeLDF();
               iewTemp.updateCren();
               iewTemp.writeCIF();
               
               % Call InductEx with necessary switches
               [status, result] = dos(['chdir ', iew.foldRoot, ' & ', ...
                  'inductex washer.cif -l UCBerkeley.ldf -i ', ...
                  'x.inp -fh -n washer.cir -k']);
               % Generated files:
               %     a.txt          = blank file
               %     fastout.out    = Current output file (InductEx)
               %     ix.cur         = Por current output file (FastHenry)
               %     j_P1.mat       = Current distribution (FastHenry)
               %     sol.txt        = Inductance + runtime (InductEx)
               %                       Also, it's the same as 'result'
               %     Zc.mat         = Complex impedance (FastHenry)
               
               str = regexpi(result,['(?<=Inductor\s*Design\s*',...
                  'Extracted\s*AbsDiff\s*PercDiff\s*)',...
                  '(\S*[ \t]*)*'],'match');
               iewTemp.L = cell2mat(textscan(str{1},'%*s %*f %f %*f %*f'));
               fprintf(1,'Extracted inductance: %g pH\n',iewTemp.L);
               
               if status
                  disp(result);
                  error('iewasher:run','Error running InductEx')
               end
               
               assert(exist([iewTemp.foldRoot,'x.inp'],'file') == 2, ...
                  'The file x.inp does not exist.')
               assert(exist([iewTemp.foldRoot,'J_P1.mat'],'file') == 2, ...
                  'iewasher:run:J_P1_missing',...
                  'The file J_P1.mat does not exist.')
               
               iewTemp.importINP
               iewTemp.importJ
               
               iewTemp.cren.mode = 'auto';
               
               if debugPlot
                  iewTemp.plotCIF
                  iewTemp.plotJ
               end
               
               % Record the current distribution:
               [xJ, J, ~] = iewTemp.JindsAcrossLinewidth;
               iewTemp.iterInfo{iewTemp.runIter+1}.JacrossW = [xJ,J];
               
               % Calculate Phi^2 for surface spins:
               if strcmpi(p.Results.calcBSurf,'iter')
                  iewTemp.calcTotB;
                  iewTemp.iterInfo{iewTemp.runIter+1}.Phi2surf = ...
                     iewTemp.Phi2surf;
                  iewTemp.checkAnalyticFun;
               end
               
               % Calculate Phi^2 for edge spins:
               if strcmpi(p.Results.calcBEdge,'iter')
                  iewTemp.calcBedge;
                  iewTemp.iterInfo{iewTemp.runIter+1}.Phi2edge = ...
                     iewTemp.Phi2edge;
               end
               
               % If all the previous code executes successfully, we can
               % destroy the old copy of iew and replace it with the new
               % one. Since the particular handle of iew is probably used
               % outside this program, it's important to preserve it and
               % just copy the data rather than creating an entirely new 
               % object. The copy function that I wrote does this, but you 
               % have to pass both handles.
               copy(iewTemp,iew);
               % The temporary handle and its associated object are no 
               % longer needed.
               delete(iewTemp);
               
               % Current iteration successful, increase the counter:
               iew.runIter = iew.runIter + 1;
            catch err % Iteration failed
               iew.err = err;
               % We've hit an error, so stop iterating and print the error
               % information:
               fprintf(1,'Oh no! I ran into an error on iteration %d:\n',...
                  iew.runIter + 1);
               disp(err.identifier)
               disp(err.message)
               for i = 1:length(err.stack)
                  fprintf('\tIn <a href="matlab:opentoline(''%s'', %d, %d)">%s</a>\n', ...
                     err.stack(i).file, err.stack(i).line, 1, err.stack(i).name);
               end
               hitErr = true;
               
               % Sometimes the code will crash if files are left open:
               fid = fopen('all');
               for i=1:length(fid), fclose(fid(i)); end
               clear fid i
            end % try/catch
         end % while
         
         % Calculate Phi^2 for surface spins:
         if strcmpi(p.Results.calcBSurf,'end')
            iew.calcTotB;
            iew.iterInfo{iew.runIter}.Phi2surf = iew.Phi2surf;
            iew.checkAnalyticFun;
         end
         
         % Calculate Phi^2 for edge spins:
         if strcmpi(p.Results.calcBEdge,'end')
            iew.calcBedge;
            iew.iterInfo{iew.runIter}.Phi2edge = iew.Phi2edge;
         end
         
         if nargout == 1
            varargout{1} = iew.objData;
         end
      end % function
      
      %% importINP
      function importINP(iew)
         % Import .inp file generated by FastHenry
         
         % Path to the input file:
         filename = [iew.foldRoot,'x.inp'];
         
         % Nodes is a struct containing relevant info about nodes:
         iew.nodes = struct('ind',zeros(0,'uint16'), ...
            'x',[],'y',[],'z',[],'A',[]);
         % Edges is a struct containing relevant info about edges:
         iew.edges = struct(...
            'ind',zeros(0,'uint16'), ... % to save memory
            'N1', zeros(0,'uint16'), ...
            'N2', zeros(0,'uint16'), ...
            'l',[],'w',[],'h',[], 'dir','');
         
         fid = fopen(filename);
         while ~feof(fid)
            curr = fgetl(fid);
            if length(curr)<1, continue, end
            switch upper(curr(1))
               case '*', continue
               case 'N' % node
                  tmp = sscanf(curr,'N%d x=%f y=%f z=%f');
                  iew.nodes.ind(end+1,1) = tmp(1);
                  iew.nodes.x(end+1,1) = tmp(2);
                  iew.nodes.y(end+1,1) = tmp(3);
                  iew.nodes.z(end+1,1) = tmp(4);
               case 'E' % edge
                  tmp = sscanf(curr,'E%d N%d N%d w=%f h=%f');
                  iew.edges.ind(end+1,1) = tmp(1);
                  iew.edges.N1(end+1,1) = tmp(2);
                  iew.edges.N2(end+1,1) = tmp(3);
                  iew.edges.w(end+1,1) = tmp(4);
                  iew.edges.h(end+1,1) = tmp(5);
                  
                  % Determine direction of edge (vertical or horizontal?):
                  if abs(iew.nodes.x(iew.edges.N1(end)) - ...
                        iew.nodes.x(iew.edges.N2(end))) < 1e-6
                     iew.edges.dir(end+1,1) = 'v'; % x coords same
                     iew.edges.l(end+1,1) = ... % compute the length
                        iew.nodes.y(iew.edges.N2(end)) - ...
                        iew.nodes.y(iew.edges.N1(end));
                  elseif abs(iew.nodes.y(iew.edges.N1(end)) - ...
                        iew.nodes.y(iew.edges.N2(end))) < 1e-6
                     iew.edges.dir(end+1,1) = 'h'; % y coords same
                     iew.edges.l(end+1,1) = ... % compute the length
                        iew.nodes.x(iew.edges.N2(end)) - ...
                        iew.nodes.x(iew.edges.N1(end));
                  else
                     iew.edges.dir(end+1) = 'x'; % neither
                     warning('iewasher:importINP:nodeWithoutDir',...
                        'Node without direction.')
                  end
            end
         end
         fclose(fid);
         
         % Find and save the midpoint of the segment (useful for
         % plotting):
         iew.edges.xMid = mean([iew.nodes.x(iew.edges.N1), ...
            iew.nodes.x(iew.edges.N2)],2);
         iew.edges.yMid = mean([iew.nodes.y(iew.edges.N1), ...
            iew.nodes.y(iew.edges.N2)],2);
         
         % Set all the node areas to nan. We'll calculate them later.
         iew.nodes.A = nan(size(iew.nodes.ind));
      end
      
      %% importJ
      function importJ(iew)
         % Import the current distribution from the exported FastHenry file
         
         % Read from text file:
         filename = [iew.foldRoot,'J_P1.mat'];
         assert(exist(filename,'file') == 2, ...
            'iewasher:importJ:J_P1_missing', ...
            'The file J_P1.mat does not exist.')
         J = dlmread(filename); %#ok<*PROP>
         J(:,[3 6]) = []; % throw away imported z components
         J(:,[1,2]) = J(:,[1,2])*1e6; % convert xy units to microns
         
         if size(J,1) > 2*length(iew.nodes.ind)
            warning('iewasher:importJ:notEnoughNodes',...
               ['It looks like FastHenry automatically increased '...
               'the resolution.\nTry decreasing gapMax.']);
         end
         
         % Associate every xy position of imported J to a corresponding
         % node in the mesh, if there exists one.
         warning('off','iewasher:importJ:notEnoughNodes')
         warning('on','iewasher:importJ:JwithNoCorrespondingNode')
         JtoNode = zeros(size(J(:,1)));
         for i = 1:length(JtoNode)
            JtoNodei = iew.nodes.ind(...
               abs(iew.nodes.x - J(i,1)) < 1e-4 & ...
               abs(iew.nodes.y - J(i,2)) < 1e-4);
            if length(JtoNodei) ~= 1
               % J point is very far from any node:
               if max([min(abs(iew.nodes.x - J(i,1))), ...
                     min(abs(iew.nodes.y - J(i,2)))]) > 0.01
                  warning('iewasher:importJ:notEnoughNodes',...
                     ['It looks like FastHenry automatically increased '...
                     'the resolution.\nTry decreasing gapMax.'])
                  warning('off','iewasher:importJ:notEnoughNodes')
               end
               warning('iewasher:importJ:JwithNoCorrespondingNode',...
                  'Current segment found without corresponding node.')
               warning('off','iewasher:importJ:JwithNoCorrespondingNode')
               continue
            end
            JtoNode(i) = JtoNodei;
         end
         
         % Throw away imported J elements that don't correspond to a
         % node in the mesh:
         J = J(JtoNode ~= 0,:);
         JtoNode = JtoNode(JtoNode ~= 0,:);
         
         % If things are working, there should be a current
         % corresponding to every segment:
         assert(size(J,1) == length(iew.edges.ind), ...
            'Not every segment has a current.');
         
         % Match imported current to the correct edge:
         iew.edges.J = zeros(length(iew.edges.ind),2);
         JtoEdge = zeros(size(J,1),1);
         for i = 1:length(JtoEdge)
            % Edge indices where the J-node matches N1:
            indi = JtoNode(i) == iew.edges.N1;
            if ~any(J(i,3:4)) % no current (can't deduce direction)
            elseif J(i,3)~=0 && J(i,4) == 0 % horizontal segment
               indi = indi & iew.edges.dir == 'h';
               assert(nnz(indi) == 1, ...
                  'Nonunity number of horizontal segments starting at N1')
               JtoEdge(i) = find(indi);
            elseif J(i,3)==0 && J(i,4) ~= 0 % vertical segment
               indi = indi & iew.edges.dir == 'v';
               assert(nnz(indi) == 1, ...
                  'Nonunity number of vertical segments starting at N1')
               JtoEdge(i) = find(indi);
            else
               warning('iewasher:importJ:noEdgeFound',...
                  'Current segment not matched with edge.')
            end
         end
         
         % Add the imported current to the edges:
         iew.edges.J(find(JtoEdge),:) = J(find(JtoEdge),3:4); %#ok<FNDSB>
         
         % For instance, we can plot the current segments, but
         % associated with the second node (N2) instead of the first:
         % figure, quiver(...
         %     iew.nodes.x(iew.edges.N2), iew.nodes.y(iew.edges.N2), ...
         %     iew.edges.J(:,1), iew.edges.J(:,2))
         
         % Sum all the currents going into each node to get node current:
         iew.nodes.J = zeros(length(iew.nodes.ind),2);
         for i = 1:length(iew.nodes.J)
            inds = iew.nodes.ind(i) == iew.edges.N1 | ...
               iew.nodes.ind(i) == iew.edges.N2;
            if any(any(iew.edges.J(inds,:)))
               iew.nodes.J(i,:) = sum(iew.edges.J(inds,:),1) ./ ...
                  sum(iew.edges.J(inds,:) ~= 0,1);
               % Divide by 0 introduces a NaN:
               iew.nodes.J(i,sum(iew.edges.J(inds,:) ~= 0,1) == 0) = 0;
            end
         end
         
         % Normalize the current to the total flowing around the loop:
         iew.normalizeJ;
         
      end %function
      
      %% normalizeJ
      function normalizeJ(iew)
         % Function normalizes the current such that the total current
         % flowing around the loop is unity.
         
         % Find indices that correspond to a slice across the bottom
         % side of the washer along x = 0 (approximately):
         inds = iew.edges.dir == 'h' & iew.edges.yMid < 0;
         x = iew.edges.xMid(inds);
         [~, xMinInd] = min(abs(x));
         inds = abs(iew.edges.xMid - x(xMinInd)) < 1e-3 & inds;
         assert(nnz(inds) > 0,'iewasher:normalizeJ:noInds',...
            'No vertical segments found to normalize current.')
         % figure, plot(iew.edges.yMid(inds),-iew.edges.J(inds,1),'o')
         
         % The total current is the sum of the y components (weighted by
         % the width of the segment):
         assert(all(~iew.edges.J(inds,2)),...
            'iewasher:normalizeJ:noVerticalCurrent', ...
            'There should be vertical current, but there is.');
         assert(abs(sum(iew.edges.w(inds)) - iew.W) < 1e-6, ...
            'iewasher:normalizeJ:widthsNotCorrect', ...
            'The widths of segments don''t add up to correct width');
         normJ = abs(sum(iew.edges.J(inds,1) .* ... % current density
            iew.edges.w(inds) .* ... % times the width
            iew.thickness )); % and the thickness;
         
         % Normalize to the total current:
         iew.nodes.J = iew.nodes.J / normJ;
         iew.edges.J = iew.edges.J / normJ;
      end
      
      %% Crenellation functions
      
      function updateCren(iew)
         % Updates the crenellation widths (either default or auto)
         
         switch upper(iew.cren.mode)
            case 'DEFAULT'
               % This case defines the default crenellation widths that are
               % used for the first iteration.
               
               % --- crenellations across linewdith ---
%                Nins = 8; % N crenellations on inside edge
               Nins = 8; % N crenellations on inside edge
               Nmid = 8; % N crenellations in middle of segment
               Nout = 4; % N crenellations on the outside edge
               
               iew.cren.line = zeros(Nins+Nmid+Nout,1);
               
               iew.cren.line(1:Nins) = (iew.W/100)*1e3;
               iew.cren.line(end-Nout+1:end) = (iew.W/50)*1e3;
               iew.cren.line = round(iew.cren.line/2)*2;
               
               iew.cren.line(Nins+1:Nins+Nmid) = ...
                  round((iew.W*1e3-sum(iew.cren.line))/Nmid);
               
               iew.cren.line(Nins+Nmid) = 0;
               iew.cren.line(Nins+Nmid) = ...
                  round(iew.W*1e3-sum(iew.cren.line));
               
               % --- crenellations across hole ---
               
               Nout = 6;  % N crenellations near the corners
               Nmid = 14; % N crenellations in the middle
               
               iew.cren.hole = zeros(2*Nout+Nmid,1);
               iew.cren.hole(1:Nout) = (iew.W/100)*1e3;
               iew.cren.hole(end-Nout+1:end) = (iew.W/100)*1e3;
               iew.cren.hole = round(iew.cren.hole/2)*2;
               
               iew.cren.hole(Nout+1:Nout+Nmid) = ...
                  round((iew.iDim(1)*1e3-sum(iew.cren.hole))/Nmid);
               
               iew.cren.hole(Nins+round(Nmid/2)) = 0;
               iew.cren.hole(Nins+round(Nmid/2)) = ...
                  round(iew.iDim(1)*1e3-sum(iew.cren.hole));
               
               assert(abs((sum(iew.cren.hole)*1e-3 - iew.iDim)/iew.iDim)...
                  < 1e-6, ...
                  'Initial hole crenellations don''t sum correctly.')
            case 'AUTO'
               % Automatically generate crenellations based on the current
               % distribution solved by Fasthenry:
               iew.calcAutoCren;
         end
      end

      function calcAutoCren(iew)
         % Calculates the automatic crenellation widths, defined based on
         % the current distribution from the previous iteration.
         
         % ----- Crenellations for the linewidth -----
         N = iew.cren.N(1);
         maxSegW = iew.W/8;
                  
         % --- From the cross-section, optimize the segment widths ---
         [y, J, ~] = iew.JindsAcrossLinewidth;
         iew.cren.line = iew.optimizeCrenWidths(...
            y, J, [-iew.oDim(1)/2, -iew.iDim(1)/2], N, maxSegW);         
         
         
         % ----- Crenellations for the hole -----
         N = iew.cren.N(2);
         maxSegW = iew.iDim(1)/10;
         
         % Find the indices of the segments just to the left of the hole:
         inds = iew.edges.dir == 'v' & iew.edges.xMid < -iew.iDim(1)/2;
         inds = inds & abs(iew.edges.yMid) < iew.iDim(2)/2 & ...
            abs(iew.edges.xMid - max(iew.edges.xMid(inds))) < 1e-3;
         assert(nnz(inds) > 0, 'iewasher:calcAutoCren:noHoleIndsFound')
         
         % --- From the cross-section, optimize the segment widths ---
         iew.cren.hole = iew.optimizeCrenWidths(...
            iew.edges.yMid(inds), iew.edges.J(inds,2), ...
            [-iew.iDim(2)/2, iew.iDim(2)/2], N, maxSegW);
         
         % There should be an odd number of crenellations across the hole.
         % This code attempts to enforce that condition.
         if mod(length(iew.cren.hole),2) == 0
            len = length(iew.cren.hole);
            inds = [len/2, len/2+1];
            % Basically, split the middle two segments into three:
            iew.cren.hole(inds) = ...
               round(sum(iew.cren.hole(inds))/3);
            iew.cren.hole = [iew.cren.hole(1:inds(1)); ...
               round(iew.iDim(1)*1e3 - sum(iew.cren.hole)); ...
               iew.cren.hole(inds(2):len)];
         end % if
      end % function

      function crenW = optimizeCrenWidths(iew,x,J,lims,N,maxSegW)
         % Given a current distribution, optimize the crenellation widths
         
         % Create a spline linear interpolation of the current density. I
         % tried other ways, but most methods don't extrapolate:
         % s1 = csapi(x,J);
         % s1 = csaps(x,J,0.8);
         s1 = interp1(x,J,'linear','pp');
         % figure, fnplt(s1), hold on, plot(x,J,'o')
         
         % Calculate the max density on the left and right, and the min:
         JmaxL = fnval(s1,lims(1));
         JmaxR = fnval(s1,lims(2));
         [Jmin, xMin] = fnmin(s1);
         assert(lims(1) < xMin && xMin < lims(2), ...
            'iewasher:optimizeCrenWidths:limitsNotCorrect', ...
            'Limits not correct.')
         
         % dJ is the maximum amount in current density that we want to
         % allow between adjacent segments. If we want roughly N segments,
         % then dJ corresponds to the total amount that J changes divided
         % by N. In practice, this gives more segments than N because I
         % also want to keep the width of segments in the middle to at most
         % some fraction of the linewidth (say, W/10).
         Jchange = (JmaxL - Jmin) + (JmaxR - Jmin);
         dJ = Jchange / N;
         
         % Now we go through and find where the edges of the segments
         % should be in order to comply with the dJ and W/10 constraints.
         % The vector z will store the x-coordinates of these boundaries.
         z = lims(1);
         
         % Starting from the left side, search for boundaries one by one
         % and add them to z. Stop the search when you get to xMin
         % because this is where the current starts to increase.
         Jnext = JmaxL - dJ;
         while isempty(z) || (z(end) < xMin)
            % Temporary spline interpolation just to find where the current
            % density equals Jnext:
            s2 = interp1(x,J - Jnext,'linear','pp');
            % Find the x-coordinate where the current equals Jnext:
            zi = fnzeros(s2,[lims(1), xMin]);
            if ~isempty(zi)
               % Choose zi unless it steps too far:
               z(end+1) = min([z(end) + maxSegW, min(zi(1,:))]);
               Jnext = Jnext - dJ;
            else
               % If no zero is found, use the maximum width allowed or
               % xMin, whichever is smaller.
               z(end+1) = min([z(end) + maxSegW, xMin]);
               Jnext = fnval(s1,z(end)) - dJ;
            end
         end
         
         % Starting from the right side, search for boundaries. Stop when
         % you get to xMin.
         Jnext = JmaxR-dJ;
         while z(end) >= xMin
            s2 = interp1(x,J - Jnext,'linear','pp');
            zi = fnzeros(s2,[xMin, lims(2)]);
            if ~isempty(zi)
               if z(end) - max(zi(1,:)) > maxSegW
                  z(end+1) = z(end) - maxSegW;
                  Jnext = fnval(s1,z(end)) - dJ;
               else
                  z(end+1) = max(zi(1,:));
                  Jnext = Jnext - dJ;
               end
            else
               % We've reached the middle so add the point and break the
               % loop:
               z(end+1) = z(end) - maxSegW;
               if z(end) < xMin
                  z(end) = xMin;
                  break
               end
               Jnext = fnval(s1,z(end)) - dJ;
            end
         end
         z = [sort(z(1:end-1)), lims(2)];
         
         % Remove boundaries that are really close together:
         zIndToCombine = diff(z)' < diff(lims)/50 & ...
            abs(diff(fnval(s1,z)')) < JmaxR/100;
         for i = find(zIndToCombine)'
            z(i) = mean(z([i,i+1]));
            z(i+1) = [];
         end
         
         % Convert from boundaries to crenellation widths (making them
         % multiples of 2 nm so that the midpoint doesn't need to be
         % rounded):
         crenW = round(diff(z)'*1e3/2)*2;
         % Reverse the ordering because crenellations are written right to
         % left:
         crenW = crenW(end:-1:1);
         % Adjust the size of the last crenellation to account for any
         % rounding errors:
         crenW(end) = round(diff(lims)*1e3) - sum(crenW(1:end-1));
         assert(sum(crenW) == round(diff(lims)*1e3), ...
            'Crenellations don''t add up to difference of limits.');
         
         if all(abs(x) >= iew.iDim(1)/2 & abs(x) <= iew.oDim(1)/2)
            iew.iterInfo{iew.runIter}.optimizeCWline.s1 = s1;
            iew.iterInfo{iew.runIter}.optimizeCWline.x = x;
            iew.iterInfo{iew.runIter}.optimizeCWline.z = z;
            iew.iterInfo{iew.runIter}.optimizeCWline.J = J;
         elseif all(abs(x) <= iew.iDim(1)/2)
            iew.iterInfo{iew.runIter}.optimizeCWhole.s1 = s1;
            iew.iterInfo{iew.runIter}.optimizeCWhole.x = x;
            iew.iterInfo{iew.runIter}.optimizeCWhole.z = z;
            iew.iterInfo{iew.runIter}.optimizeCWhole.J = J;
         end
         
         if iew.plotOptimizeCrenWidths
            figure
            fnplt(s1), hold on
            plot(x, J,'o')
            plot(z, fnval(s1,z),'r+')
         end
      end

      function plotOptimizedCrenWidths(iew)
         % Plot the optimized crenellation widths, along with current
         % density, for each of the iterations
         
         figure
         for i = 1:length(iew.iterInfo)
            if isfield(iew.iterInfo{i},'optimizeCWline')
               s1 = iew.iterInfo{i}.optimizeCWline.s1;
               x = iew.iterInfo{i}.optimizeCWline.x;
               z = iew.iterInfo{i}.optimizeCWline.z;
               J = iew.iterInfo{i}.optimizeCWline.J;
               fnplt(s1), hold on
               plot(x, J,'o')
               plot(z, fnval(s1,z),'r+')
            end
         end
         
         figure
         for i = 1:length(iew.iterInfo)
            if isfield(iew.iterInfo{i},'optimizeCWhole')
               s1 = iew.iterInfo{i}.optimizeCWhole.s1;
               x = iew.iterInfo{i}.optimizeCWhole.x;
               z = iew.iterInfo{i}.optimizeCWhole.z;
               J = iew.iterInfo{i}.optimizeCWhole.J;
               fnplt(s1), hold on
               plot(x, J,'o')
               plot(z, fnval(s1,z),'r+')
            end
         end
         
      end
      
      %% calcNodeAreas
      function calcNodeAreas(iew,varargin)
         % Function calculates the area corresponding to each node
         
         %{
         This function calculates the effective node areas for the nodes in
         the lower left quarter of the washer. A voronoi diagram
         essentially does this, but I didn't figure out how to do boundary
         conditions on it, so you have to manually edit all the point on
         the edges. (Is there a better way to do this?)
         %}
         
         % ----- Interpret the input flags -----
         debugPlot = false;
         i = 1;
         while i < nargin
            switch upper(varargin{i})
               case 'PLOT'
                  debugPlot = true;
                  i = i+1;
            end
         end
         
         % Indices of nodes in bottom left quarter of the washer,
         % but not in a crenellation:
         inds = iew.nodes.x <= 0 & iew.nodes.y <= 0 & ...
            iew.nodes.x > -iew.oDim(1)/2 & iew.nodes.y > -iew.oDim(2)/2;
         
         % Coordinates of nodes that we want to keep
         xn = iew.nodes.x(inds);
         yn = iew.nodes.y(inds);
         
         % Voronoi indices, which we'll use to compute the effective area
         % that each node represents:
         [v,c] = voronoin([xn, yn]);
         
         if debugPlot, iew.plotCIF, iew.plotNodes, end
         
         % Find the voronoi points that lie inside the hole and outside the
         % washer:
         vToRem = (v(:,1) > -iew.iDim(1)/2 & v(:,2) > -iew.iDim(2)/2) | ...
            (v(:,1) <= -iew.oDim(1)/2 | v(:,2) <= -iew.oDim(2)/2);
         vToRemInds = find(vToRem);
         
         % Some useful values:
         xnMin = min(xn);                ynMin = min(yn);
         xnMax = max(xn);                ynMax = max(yn);
         xnIns = max(xn(yn == max(yn)));   ynIns = max(yn(xn == max(xn)));
         
         % Function forms a square out of the (x,y) coordinates:
         cSq = @(xn,yn) [min(xn),min(yn); max(xn),min(yn); ...
            max(xn),max(yn); min(xn),max(yn)];
         
         % Loop over all nodes:
         nodeAreas = zeros(length(c),1);
         for i = 1:length(c)
            % Remove any points that may be in the hole or out at infinity:
            [~, ia, ~] = intersect(c{i}, vToRemInds);
            c{i}(ia) = [];
            
            switch length(c{i})
               case 4 % easy case: all nodes specified by voronoi points:
                  vtmp = v(c{i},:);
               case 2 % along an edge somewhere
                  if xn(i) == xnMin % node is on left side
                     vtmp = cSq([v(c{i},1);-iew.oDim(1)/2], v(c{i},2));
                  elseif xn(i) == xnIns && yn(i) >= ynIns % left inside
                     vtmp = cSq([v(c{i},1);-iew.iDim(1)/2], v(c{i},2));
                  elseif xn(i) == xnMax % right side
                     vtmp = cSq([v(c{i},1);             0], v(c{i},2));
                     
                  elseif yn(i) == ynMin % node is on bottom
                     vtmp = cSq(v(c{i},1), [v(c{i},2);-iew.oDim(2)/2]);
                  elseif yn(i) == ynIns && xn(i) >= xnIns % bottom inside
                     vtmp = cSq(v(c{i},1), [v(c{i},2);-iew.iDim(2)/2]);
                  elseif yn(i) == ynMax % top side
                     vtmp = cSq(v(c{i},1), [v(c{i},2);             0]);
                  else
                     error('iewasher:calcTotB','Node not recognized.')
                  end % if
                  
               case 3 % two nodes in the inside corner
                  if max(v(c{i},1)) < -iew.iDim(1)/2 && yn(i) > ...
                        -iew.iDim(2)/2
                     vtmp = [v(c{i},:); -iew.iDim/2; ...
                        -iew.iDim(1)/2, max(v(c{i},2))];
                     vtmp = vtmp(convhull(vtmp),:);
                     vtmp = vtmp(1:end-1,:);
                  elseif max(v(c{i},2)) < -iew.iDim(2)/2 && xn(i) > ...
                        -iew.iDim(1)/2
                     vtmp = [v(c{i},:); -iew.iDim/2; ...
                        max(v(c{i},1)), -iew.iDim(2)/2];
                     vtmp = vtmp(convhull(vtmp),:);
                     vtmp = vtmp(1:end-1,:);
                  else
                     vtmp = cSq(v(c{i},1), v(c{i},2));
                     nInHole = find(...
                        -iew.iDim(1)/2 < vtmp(:,1) - 1e-6 & ...
                        -iew.iDim(2)/2 < vtmp(:,2) - 1e-6);
                     if nInHole == 3
                        vtmp = [...
                           [1 0 0 0 0 1]' * min(vtmp(:,1)) + ...
                           [0 1 1 0 0 0]' * max(vtmp(:,1)) + ...
                           [0 0 0 1 1 0]' * -iew.iDim(1)/2, ...
                           [1 1 0 0 0 0]' * min(vtmp(:,2)) + ...
                           [0 0 0 0 1 1]' * max(vtmp(:,2)) + ...
                           [0 0 1 1 0 0]' * -iew.iDim(2)/2];
                     end
                     
                  end
                  
               case 1 % nodes on outside corners
                  if xn(i) == xnMin && yn(i) == ynMin
                     vtmp = cSq([v(c{i},1);-iew.oDim(1)/2], ...
                        [v(c{i},2);-iew.oDim(1)/2]);
                  elseif xn(i) == xnMin && yn(i) == ynMax
                     vtmp = cSq([v(c{i},1);-iew.oDim(1)/2], ...
                        [v(c{i},2);0]);
                  elseif xn(i) == xnIns && yn(i) == ynMax
                     vtmp = cSq([v(c{i},1),-iew.iDim(1)/2], [v(c{i},2),0]);
                  elseif xn(i) == xnMax && yn(i) == ynMin
                     vtmp = cSq([v(c{i},1),0], [v(c{i},2),-iew.oDim(2)/2]);
                  elseif xn(i) == xnMax && yn(i) == ynIns
                     vtmp = cSq([v(c{i},1),0], [v(c{i},2),-iew.iDim(2)/2]);
                  else
                     error('iewasher:calcTotB','Node not recognized.')
                  end
            end % switch
            
            K = [1:size(vtmp,1), 1];
            nodeAreas(i) = polyarea(vtmp(K,1),vtmp(K,2));
            if debugPlot, plot(vtmp(K,1),vtmp(K,2)); end
         end % for
         
         if any(nodeAreas == 0)
            warning('iewasher:calcTotB:nodeAreaEqualsZero', ...
               'At least one node area equals zero.')
         end
         
         % nodeAreas adds up to correct total area:
         correctTotArea = (prod(iew.oDim)-prod(iew.iDim))/4;
         relAreaError = abs(correctTotArea-sum(nodeAreas))/correctTotArea;
         assert(relAreaError < 1e-4,'iewasher:calcTotB',...
            'Node areas do not add correctly. Relative error is: %g',...
            relAreaError)
         if relAreaError > 1e-7
            warning('iewasher:calcTotB:incorrectNodeArea', ...
               'Node areas do not add up correctly. Relative error is %g',...
               relAreaError);
         end
         
         iew.nodes.A(inds) = nodeAreas;
      end
      
      %% Magnetic field functions
      
      function B = calcB(iew,x,y,z)
         % Calculate the magnetic field at arbitrary position
         
         %{
         This function calculates the magnetic field at an arbitrary
         position by summing the contribution of each current segment in
         the mesh.
         %}
         x = reshape(x,numel(x),1);
         y = reshape(y,numel(y),1);
         z = reshape(z,numel(z),1);
         
         % Pull out the relevant quantities (makes code a little easier to
         % read; may also execute faster because I'm not sure if it takes a
         % long time to access data from structs):
         
         % Centers of cubes:
         x0 = iew.edges.xMid;
         y0 = iew.edges.yMid;
         
         % Cube x, y, and z double widths:
         a = iew.edges.w/2;
         b = (iew.nodes.y(iew.edges.N2) - iew.nodes.y(iew.edges.N1)) / 2;
         c = iew.thickness/2;
         
         % Current density at each cube:
         J = iew.edges.J;
         
         % Initialize the integral K to zero:
         K = 0;
         % Loop over all vertical segments and add their contribution to
         % the integral:
         for i = find(iew.edges.dir == 'v')'
            K = K + J(i,2) * ...
               KfromCube(x-x0(i), y-y0(i), z, a(i), b(i), c);
         end
         % To get the magnetic field, we have to cross the current density
         % into the result of the integral. We could wait until the end to
         % do this because all the current densities point in the same
         % direction, that is, in the yhat direction.
         B = cross(repmat([0 1 0],size(K,1),1), K);
         
         % Now we do the same thing with the horizontal segments, but first
         % we have to recompute the size of the cubes:
         a = (iew.nodes.x(iew.edges.N2) - iew.nodes.x(iew.edges.N1)) / 2;
         b = iew.edges.w/2;
         K = 0;
         for i = find(iew.edges.dir == 'h')'
            K = K + J(i,1) * ...
               KfromCube(x-x0(i), y-y0(i), z, a(i), b(i), c);
         end
         
         % Final magnetic field is the sum of the field due to the
         % horizontal and the field due to the vertical components:
         B = B + cross(repmat([1 0 0],size(K,1),1), K);
         
         % Convert to Tesla:
         % Factor of 0.1 comes from doing
         % (mu0/4pi) * iiint[ dV (J(r) x dr)/r^3 ] --> 0.1 Tesla
         % In other words, multiply the result by 0.1 to put B in units of
         % tesla.
         B = 0.1*B;
      end
      
      function B = calcBpar(iew,x,y,z)
         % Calculate the magnetic field at arbitrary position in parallel
         
         %{
          This function calculates the magnetic field at an arbitrary
          position by summing the contribution of each current segment in
          the mesh.
          Parallel version:
          This version distributes the calculation of the magnetic field
          between however many workers are established (via matlabpool).
         %}
         
         % Initialize matlabpool if not open yet:
         if matlabpool('size') == 0
            matlabpool
         end
         
         x = reshape(x,numel(x),1);
         y = reshape(y,numel(y),1);
         z = reshape(z,numel(z),1);
         numB = length(x.*y.*z);
         
         % Pull out the relevant quantities (makes code a little easier to
         % read; may also execute faster because I'm not sure if it takes a
         % long time to access data from structs):
         
         % Centers of cubes:
         x0 = iew.edges.xMid;
         y0 = iew.edges.yMid;
         
         % Current density at each cube:
         J = iew.edges.J;
         
         % Start segment of code that is divided among workers:
         spmd
            % Loop over all vertical segments and add their contribution to
            % the integral:
            
            % Cube x, y, and z double widths:
            a = iew.edges.w/2;
            b = (iew.nodes.y(iew.edges.N2) - iew.nodes.y(iew.edges.N1)) / 2;
            c = iew.thickness/2;
            
            % Initialize the integral K to zero:
            Ki = zeros(numB,3);
            for i = drange(find(iew.edges.dir == 'v')')
               Ki = Ki + J(i,2) * ...
                  KfromCube(x-x0(i), y-y0(i), z, a(i), b(i), c);
            end
            % To get the magnetic field, we have to cross the current density
            % into the result of the integral. We could wait until the end to
            % do this because all the current densities point in the same
            % direction, that is, in the yhat direction.
            Bi = cross(repmat([0 1 0],size(Ki,1),1), Ki);
            
            % Now we do the same thing with the horizontal segments, but first
            % we have to recompute the size of the cubes:
            a = (iew.nodes.x(iew.edges.N2) - iew.nodes.x(iew.edges.N1)) / 2;
            b = iew.edges.w/2;
            c = iew.thickness/2;
            
            Ki = zeros(numB,3);
            for i = drange(find(iew.edges.dir == 'h')')
               Ki = Ki + J(i,1) * ...
                  KfromCube(x-x0(i), y-y0(i), z, a(i), b(i), c);
            end
            
            % Final magnetic field is the sum of the field due to the
            % horizontal and the field due to the vertical components:
            Bi = Bi + cross(repmat([1 0 0],size(Ki,1),1), Ki);
         end % spmd
         % Bi is a composite variable that contains the summed magnetic
         % field that each worker computed. To get the total field, we
         % have to sum the contributions that each worker calculated:
         B = zeros(length(x.*y.*z),3);
         for i=1:length(Bi), B = B + Bi{i}; end
         
         % Convert to Tesla:
         % Factor of 0.1 comes from doing
         % (mu0/4pi) * iiint[ dV (J(r) x dr)/r^3 ] --> 0.1 Tesla
         % In other words, multiply the result by 0.1 to put B in units of
         % tesla.
         B = 0.1*B;
      end
      
      function calcTotB(iew)
         % Calculate the average magnetic field and multiply by the total
         % area
         
         %{
         This function calculates the average magnetic field over the
         surface of the washer. To do this, we evaluate the field at the
         (x,y) locations of each node, which preserves the nonuniform
         meshing that we so desire. To save computation time, we do this
         only in the lower-left-hand portion and multiply by 4. Finally, we
         sum the magnetic field at each node, weighted by the node area.
         
         NOTE: This function actually returns <B^2>*A, where A is the total
         area of the washer.
         %}
         
         iew.calcNodeAreas();
         
         % Coordinates of nodes that we want to keep:
         inds = ~isnan(iew.nodes.A);
         nodeAreas = iew.nodes.A(inds);
         xn = iew.nodes.x(inds);
         yn = iew.nodes.y(inds);
         
         % --------------- Calculate the surface field ---------------
         
         z = iew.thickness/2;
         % The fudge factor is needed for numerical stability. Under
         % certain circumstances, some terms in the analytic expression
         % actually cancel analytically, but evaluate to nan or inf
         % numerically. Adding an inconsequential nudge seems to help, but
         % it doesn't really change the final answer.
         fudgeFact = 1e-6*iew.oDim(1)/2;
         
         if iew.runParallel
            B = iew.calcBpar(xn+fudgeFact, yn+fudgeFact, z+fudgeFact);
         else
            B = iew.calcB(xn+fudgeFact, yn+fudgeFact, z+fudgeFact);
         end
         
         % quiver3(xn, yn, z*ones(size(xn)), B(:,1), B(:,2), B(:,3), 'k',...
         %    'autoscalefactor',0.1)
         
         % B is in units of Tesla:
         iew.Bsurf = [xn, yn, z*ones(size(xn)), B];
         % B2 is in units of Tesla^2*um^2. Factor of 4 because nodeAreas
         % only covers one quarter of the loop. Another factor of 2 because
         % there are spins on top and bottom.
         iew.B2surf = 4 * 2 * [ ...
            sum(nodeAreas.*sum(B(:,1).^2,2)), ...
            sum(nodeAreas.*sum(B(:,2).^2,2)), ...
            sum(nodeAreas.*sum(B(:,3).^2,2)) ];
      end
      
      function calcBedge(iew)
         % Calculate contribution of edge spins
         
         %{
         This function calculates the contribution to flux noise due to
         spins on the edge of the film. That is, the spins spaced
         vertically along the thickness of the film. It does so by
         computing the magnetic field at Nz points distributed uniformly
         across the width of the film, but distributed non-uniformly along
         the length, at the same locations as the midpoints of segments.
         %}
         
         Nz = 20;
         %% --------------- Inside edge ---------------
         % Indices of edges in the lower quarter:
         inds = abs(iew.edges.yMid) < iew.iDim(2)/2 & iew.edges.xMid < 0;
         % x-value of segments bordering the hole:
         xBorderingHole = max(iew.edges.xMid(inds));
         
         inds = iew.edges.xMid == xBorderingHole & ...
            iew.edges.yMid < 0 & ...
            iew.edges.yMid > -iew.iDim(2)/2;
         
         % Coordinates where we want to evaluate the magnetic field; here,
         % I've added a point at y = 0, which we need to keep track of.
         % Since the coordinates are nonuniformly spaced, we need to
         % compute the area of each "patch" that the node represents.
         x = -iew.iDim(1)/2;
         y = [iew.edges.yMid(inds); 0];
         z = linspace(0,1,Nz)*iew.thickness/2;
         
         % y-length of all the segments (it's ugly, but it works, I think):
         yLenL = diff([-iew.iDim(2)/2;y])/2;
         yLenL(1) = 2*yLenL(1);
         yLenR = [diff(y)/2; 0];
         yLen = yLenL + yLenR;
         assert(abs(sum(yLen)-iew.iDim(2)/2)/iew.iDim(2) < 1e-6, ...
            'Length of y-segments for edge-spin calculation is incorrect.')

         zLen = [0,diff(z)/2]' + [diff(z)/2,0]';
         
         nodeAreas = yLen*zLen';
         
         [X, Y, Z] = meshgrid(x, y, z);
         xn = X(:); yn = Y(:); zn = Z(:);
         
         fudgeFact = 1e-6*iew.oDim(1)/2;
         if iew.runParallel
            B = iew.calcBpar(xn+fudgeFact, yn+fudgeFact, zn+fudgeFact);
         else
            B = iew.calcB(xn+fudgeFact, yn+fudgeFact, zn+fudgeFact);
         end
         
         % From symmetry, you get a factor of 2 by reflecting across z = 0,
         % and a factor of 4 for each of the 4 sides.
         iew.B2insEdge = 2 * 4 * sum(repmat(nodeAreas(:),1,3) .* B.^2);
         % You also get a factor of 2 by reflecting across y = x, but you
         % have to reverse the x and y components. We could calculate, but
         % I'm lazy.
         iew.B2insEdge = iew.B2insEdge(:,1:3) + iew.B2insEdge(:,[2 1 3]);
         
         %%    Outside edge
         inds = abs(iew.edges.yMid) < iew.iDim(2)/2 & ...
            iew.edges.xMid < 0 & iew.edges.xMid > -iew.oDim(1)/2;
         xBorderingOutside = min(iew.edges.xMid(inds));
         
         inds = iew.edges.xMid == xBorderingOutside & ...
            iew.edges.yMid < 0;
         
         % Coordinates where we want to evaluate the magnetic field; here,
         % I've added a point at y = 0, which we need to keep track of.
         % Since the coordinates are nonuniformly spaced, we need to
         % compute the area of each "patch" that the node represents.
         x = -iew.oDim(1)/2;
         y = [iew.edges.yMid(inds); 0];
         z = linspace(0,1,Nz)*iew.thickness/2;
         
         % y-length of all the segments (it's ugly, but it works, I think):
         yLenL = diff([-iew.oDim(2)/2;y])/2;
         yLenL(1) = 2*yLenL(1);
         yLenR = [diff(y)/2; 0];
         yLen = yLenL + yLenR;
         assert(abs(sum(yLen)-iew.oDim(2)/2)/iew.oDim(2) < 1e-6, ...
            'Length of y-segments for edge-spin calculation is incorrect.')
         % Width of segments in z-direction:
         zLen = [0,diff(z)/2]' + [diff(z)/2,0]';
         % Area of nodes:
         nodeAreas = yLen*zLen';
         
         [X, Y, Z] = meshgrid(x, y, z);
         xn = X(:); yn = Y(:); zn = Z(:);
         
         fudgeFact = 1e-6*iew.oDim(1)/2;
         if iew.runParallel
            B = iew.calcBpar(xn+fudgeFact, yn+fudgeFact, zn+fudgeFact);
         else
            B = iew.calcB(xn+fudgeFact, yn+fudgeFact, zn+fudgeFact);
         end
         
         % From symmetry, you get a factor of 2 by reflecting across z = 0,
         % and a factor of 4 for each of the 4 sides.
         iew.B2outEdge = 2 * 4 * sum(repmat(nodeAreas(:),1,3) .* B.^2);
         % You also get a factor of 2 by reflecting across y = x, but you
         % have to reverse the x and y components. We could calculate, but
         % I'm lazy.
         iew.B2outEdge = iew.B2outEdge(:,1:3) + iew.B2outEdge(:,[2 1 3]);
      end
      
      function visualizeB(iew,visType)
         % Functions to visualize the magnetic field ('surface-field' or
         % 'cross-section')
         
         switch upper(visType)
            case 'SURFACE-FIELD'
               % ----- Surface field -----
               % This code calculates the magnetic field at the top surface
               % of the washer and plots it.
               
               % Points at which the field is to be evaluated:
               x = linspace(-1.1*iew.oDim(1)/2, 1.1*iew.oDim(1)/2, 30);
               y = linspace(-1.1*iew.oDim(2)/2, 1.1*iew.oDim(2)/2, 30);
               z = iew.thickness/2;
               [X,Y] = meshgrid(x,y);
               Z = z*ones(size(X));
               
               % Calculate the field at those points:
               B = iew.calcBpar(X(:),Y(:),z);
               % Plot the washer, then plot the magnetic field vectors:
               iew.plotCIF;
               quiver3(X(:), Y(:), Z(:), B(:,1), B(:,2), B(:,3), 'b')
               
            case 'CROSS-SECTION'
               % ----- Cross section -----
               % This code calculates the magnetic field in a vertical
               % slice that goes through the washer film.
               
               % Points at which the magnetic field is to be evaluated:
               x = linspace(-1.1*iew.oDim(1)/2, -0.9*iew.iDim(1)/2, 30);
               y = 0;
               z = linspace(-4*iew.thickness/2, 4*iew.thickness/2, 30);
               [X,Z] = meshgrid(x,z);
               Y = y*ones(size(X));
               
               % Calculate the magnetic field at specified points:
               B = iew.calcBpar(X(:),y,Z(:));
               
               doPlotCIF = false;
               if doPlotCIF
                  % Calculate the field, then plot it:
                  iew.plotCIF;
                  quiver3(X(:), Y(:), Z(:), B(:,1), B(:,2), B(:,3), 'b')
                  xlabel('x ({\mu}m)')
                  ylabel('y ({\mu}m)')
                  zlabel('z ({\mu}m)')
               end;
               
               % Plot the cross section of the film:
               patch(-[1 0 0 1 1]*iew.R-[0 1 1 0 0]*iew.iDim(1)/2,...
                  [-1 -1 1 1 -1]*iew.thickness/2, ([0,176,240]+1)/256)
               hold on
               % Plot the vector field:
               quiver(X(:), Z(:), B(:,1), B(:,3), 'b')
               % axis equal
               xlabel('x ({\mu}m)'), ylabel('z ({\mu}m)')
         end % switch
      end
      
      %% Functions related to analytic answer
      function analyticAns = analyticFun(iew,varargin)
         % Calculate the analytic prediction for <Phi^2>
         
         p = inputParser;
         p.CaseSensitive = false;
         p.addParamValue('a',0.27,@isnumeric);
         p.addParamValue('R','out', ...
            @(x) any(strcmpi(x,{'ins','out','mid'})) );
         p.parse(varargin{:});
         
         a = p.Results.a;
         Rstr = ['R',p.Results.R];
         
         % To vectorize and return in the same shape as iew:
         R = reshape([iew.(Rstr)],size(iew));
         W = reshape([iew.W],size(iew));
         thickness = reshape([iew.thickness],size(iew));
         lambda = reshape([iew.lambda],size(iew));
         
         % ----- Analytic formula from Bialczak et al. -----
         % Wolfram code to get the correct prefactor:
         % (2*mu0^2/3) * (bohr magneton)^2 * (5E17/m^2) / (flux quantum)^2
         analyticAns = 1.058772e-11 * R./W .* ...
            (log(2*thickness.*W./lambda.^2)/(2*pi) + a);
         
         % "Convert" from a circle to a square geometry:
         analyticAns = analyticAns * 4/pi;
      end
      
      function J = analyticJacrossW(iew,x,varargin)
         % Calculate normalized analytic current distribution across
         % linewidth (J/um^2)
         
         % Interpret input arguments:
         p = inputParser;
         p.CaseSensitive = false;
         p.addParamValue('a',exp(-0.27*2*pi),@isnumeric); % from Bialczak
         p.parse(varargin{:});
         
         % Define relevant variables:
         a = p.Results.a;
         w = iew.W;
         lam = iew.lambda;
         b = iew.thickness;
         
         % Define the function and crossover point x0
         x0 = w/2 - a*lam^2/(2*b);
         function val = analyticJ(w,b,lam,a,x0,x)
            v1 = (1-(2*x/w).^2).^(-1/2) .* (abs(x) <= x0);
            v1(~isfinite(v1)) = 0;
            
            v2norm = (1-(2*x0/w).^2).^(-1/2) / ...
               exp(-(w/2-abs(x0))*b/(a*lam^2));
            % Definition given in Van Duzer that leads to discontinuity:
            % v2norm = (1.165/lam) * (w*b/a)^(1/2);
            v2 = v2norm * exp(-(w/2-abs(x))*b/(a*lam^2)) .* (abs(x) > x0);
            val = v1+v2;
         end
         
         % Compute the current density (A/um)
         J = analyticJ(w,b,lam,a,x0,x) / ...
            quad(@(xx) analyticJ(w,b,lam,a,x0,xx),-x0,x0);
         % Convert to A/um^2:
         J = J / iew.thickness;
      end
      
      function checkAnalyticFun(iew,varargin)
         % Check the accuracy of the analytic answer by comparing it to the
         % numerical answer.
         
         % Parse the input arguments
         p = inputParser;
         p.CaseSensitive = false;
         p.addOptional('plot',false);
         p.addParamValue('a',exp(-0.27*2*pi),@isnumeric); % from Bialczak
         p.parse(varargin{:});
         
         plotOn = p.Results.plot;
         a = p.Results.a;
         
         % ----- Compare current distributions -----
         [yJ, J, inds] = iew.JindsAcrossLinewidth;
         
         w = iew.W;
         % Analytic function as givin in Van Duzer:         
         function val = JfunFull(w,b,lam,a,x)
            x0inFun = w/2 - a*lam^2/(2*b);
            v1 = (1-(2*x/w).^2).^(-1/2) .* (abs(x) <= x0inFun);
            v1(~isfinite(v1)) = 0;
            
            v2norm = (1-(2*x0inFun/w).^2).^(-1/2) / ...
               exp(-(w/2-abs(x0inFun))*b/(a*lam^2));
            % v2norm = (1.165/lam) * (w*b/a)^(1/2);
            v2 = v2norm * exp(-(w/2-abs(x))*b/(a*lam^2)) .* ...
               (abs(x) > x0inFun);
            val = v1+v2;
         end
         
         x0 = w/2 - a*iew.lambda^2/(2*iew.thickness);
         Jfun = @(x) JfunFull(w, iew.thickness, iew.lambda, a, x);
         
%          JnormFactor = 1 / (quad(Jfun, -w/2, w/2) * iew.thickness);
         JnormFactor = 1 / (quad(Jfun, -x0, x0) * iew.thickness);
         
         warning('off','MATLAB:rankDeficientMatrix');
         coeff = nlinfit(yJ+iew.R-w/2, J, ...
            @(c,x) c(1)*JfunFull(w,iew.thickness,iew.lambda,c(2),x), ...
            [JnormFactor, 0.5]);
         iew.aOptimum = coeff(2);
         warning('on','MATLAB:rankDeficientMatrix');
         
         % Units of A^2/um^3
         % Janalytic = JnormFactor*Jfun(yJ+iew.R-w/2);
%          J2analytic = quad(@(x) (JnormFactor*Jfun(x)).^2, -w/2, w/2);
         J2analytic = quad(@(x) (JnormFactor*Jfun(x)).^2, -x0, x0);
         J2numeric = sum(J.^2.*iew.edges.w(inds));
         iew.relJ2AnalyticNumeric = J2analytic/J2numeric;
      
         J2b = J2analytic * iew.thickness^2; % Units of A^2/um
         % quad(@(x) (JnormFactor*Jfun(x)*iew.thickness).^2, -w/2, w/2);
         % Wolfram code for conversion factor for int(J^2) to <Phi^2>
         % (based on Eq. (5) from Bialczak)
         % (4/pi)*(pi/6)*(mu_0)^2*(bohr magneton)^2*(5E17/m^2)*1um*(1/um)/(flux quantum)^2
         % Factor of 1um comes from R
         % Factor if (1/um) comes from int(J^2)/int(J)^2
         J2bToPhi2 = 1.058772e-11;
         iew.Phi2estimateFromAnalyticJ = J2bToPhi2 * J2b * iew.R;
                  
         if plotOn
            figure, subplot(2,2,1)
            x = linspace(-w/2,w/2,10000);
            plot(yJ,J), hold on, plot(x-iew.R+w/2,JnormFactor*Jfun(x),'r')
            plot(x-iew.R+w/2,JnormFactor*JfunFull(w, iew.thickness, iew.lambda, iew.aOptimum, x),'k')
            xlabel('x ({\mu}m)')
            ylabel('J [(A/({\mu}m)^2]')
            legend({'Numerical','Analytic','Analytic, aOptimum'})
            
            subplot(2,2,2)
            semilogx(-(yJ+iew.R-iew.W),J), hold on
            semilogx(x+w/2,JnormFactor*Jfun(x),'r')
            semilogx(x+w/2,JnormFactor*JfunFull(w, iew.thickness, iew.lambda, iew.aOptimum, x),'k')
            xlim([min(x(x>0)), w/2])
            xlabel('R-x ({\mu}m)')
            ylabel('J [(A/({\mu}m)^2]')
         end
         
         % ----- Compare magnetic fields -----
         [yB, B, ~] = iew.BindsAcrossLinewidth;
         
         % Ampere's rule to convert from (current density)*(thickness) to T:
         % (Wolfram code: (mu_0/2) * (1A/um^2) * (1um)) = 0.6283 Tesla
         JtoB = pi/5;
         
         if plotOn
            subplot(2,2,3)
            plot(yB, B(:,2), '-'), hold on
            plot(x-iew.R+w/2, JtoB*JnormFactor*Jfun(x)*iew.thickness, 'r')
            plot(x-iew.R+w/2, JtoB*JnormFactor*JfunFull(w, iew.thickness, iew.lambda, iew.aOptimum, x)*iew.thickness, 'k')
            xlabel('x ({\mu}m)')
            ylabel('B (T)')
            title('Comparison of numerical and analytic magnetic fields')
            legend({'Numerical','Analytic','Analytic, aOptimum'})
            
            subplot(2,2,4)
            semilogx(-(yB+iew.R-iew.W), B(:,2), '-'), hold on
            semilogx(x+w/2, JtoB*JnormFactor*Jfun(x)*iew.thickness, 'r')
            semilogx(x+w/2, JtoB*JnormFactor*JfunFull(w, iew.thickness, iew.lambda, iew.aOptimum, x)*iew.thickness, 'k')
            xlabel('R-x ({\mu}m)')
            ylabel('B (T)')
            xlim([min(x(x>0)), w/2])
         end
         
         iew.errorFromEpsExpansion =  ...
            2*log(2*iew.thickness*w/(a*iew.lambda^2))/(pi^2*w) / ...
            (atanh(2*x0/w)/(w*asin(2*x0/w)^2)) - 1;
         
         B2totAnalytic = quad(@(x) ...
            (JtoB*JnormFactor*Jfun(x)*iew.thickness).^2, -w/2, w/2);
         B2totAnalyticTox0 = quad(@(x) ...
            (JtoB*JnormFactor*Jfun(x)*iew.thickness).^2, -x0, x0);
         
         % Conversion factor:
         % ((5E17/m^2)*(bohr magneton)^2/3A^2)*(1T^2*um)*(8*1um)/flux quantum^2
         iew.Phi2estimateFromAnalyticB = 2.6819e-11*B2totAnalyticTox0*iew.R;
         
         B2totNumeric = sum(B(:,2).^2 .* iew.edges.w(inds));
         %   iew.relB2AnalyticNumeric = B2totNumeric / B2totAnalytic;
         iew.relB2AnalyticNumeric = B2totAnalyticTox0 / B2totNumeric;
      end
      
      %% find indices
      function [y, J, inds] = JindsAcrossLinewidth(iew,varargin)
         % Look up the current density across the linewidth
         
         p = inputParser;
         p.CaseSensitive = false;
         p.addOptional('plot',false);
         p.parse(varargin{:});
         plotOn = p.Results.plot;
         
         % Find indices that correspond to a slice across the bottom
         % side of the washer along x = 0 (approximately):
         inds = iew.edges.dir == 'h' & iew.edges.yMid < 0;
         x = iew.edges.xMid(inds);
         [~, xMinInd] = min(abs(x));
         inds = abs(iew.edges.xMid - x(xMinInd)) < 1e-3 & inds;
         assert(nnz(inds) > 0,'iewasher:normalizeJ:noInds',...
            'No vertical segments found to normalize current.')
         y = iew.edges.yMid(inds);
         J = -iew.edges.J(inds,1);
         
         if plotOn
            figure, plot(y,J,'-o')
            xlabel('x ({\mu}m)')
            ylabel('J [(A/({\mu}m)^2]')
            title('Current density across linewidth')
         end
      end
      
      function [y, B, inds] = BindsAcrossLinewidth(iew,varargin)
         % Look up the surface magnetic field across the linewidth 
         
         p = inputParser;
         p.CaseSensitive = false;
         p.addOptional('plot',false);
         p.parse(varargin{:});
         plotOn = p.Results.plot;
         
         % Find indices that correspond to a slice across the bottom
         % side of the washer along x = 0 (approximately):
         inds = iew.Bsurf(:,2) < 0;
         x = iew.Bsurf(inds,1);
         [~, xMinInd] = min(abs(x));
         inds = abs(iew.Bsurf(:,1) - x(xMinInd)) < 1e-3 & inds;
         y = iew.Bsurf(inds,2);
         B = iew.Bsurf(inds,4:6);
         
         if plotOn
            figure, plot(y, B(:,2), '-o')
            xlabel('x ({\mu}m)')
            ylabel('B (T)')
            title('Magnetic field across linewidth')
         end
      end
      
      %% Plot the current distribution
      
      function plotJ(iew,varargin)
         % Plot the current distribution at each node
         quiver(iew.nodes.x, iew.nodes.y, ...
            iew.nodes.J(:,1), iew.nodes.J(:,2), 'autoscalefactor',0.1, ...
            varargin{:})
         hold on
      end
      
      function plotJcomp(iew,varargin)
         % Plot the current distribution components. That is, plot x and y
         % components separately, just as they are solved:
         quiver(iew.edges.xMid, iew.edges.yMid, ...
            iew.edges.J(:,1), iew.edges.J(:,2),varargin{:})
         hold on
      end
      
      function plotJacrossW(iew)
         % Plot the current distribution across the linewidth as a function
         % of iteration:
      
         % Choose colors of the lines
         colors = copper(length(iew.iterInfo));
         colors = colors(end:-1:1,:);
         legendstr = {};
         
         figure
         % Loop through all the iterations and plot the current:
         for i = 1:length(iew.iterInfo)
            plot(iew.iterInfo{i}.JacrossW(:,1), ...
               iew.iterInfo{i}.JacrossW(:,2), ...
               'color',colors(i,:))
            hold on
            legendstr{i} = sprintf('Iteration %d',i);
         end
         xlabel('x ({\mu}m)')
         ylabel('J [(A/({\mu}m)^2]')
         title('Current distribution versus iteration number')
         
         legend(legendstr,'location','best');
      end
      
      function plotJbarsAcrossW(iew,varargin)
         % Plots the current distribution across W using bars:
         
         % Find info for the current distribution:
         [x, J, inds] = iew.JindsAcrossLinewidth;
         w = iew.edges.w(inds);
         
         % Plot each component:
         for i=1:length(x)
            if nargin > 1
               ColorSpecs = varargin;
            else
               ColorSpecs = {'b'};
            end
            fill(x(i) + w(i)/2*[1 1 -1 -1 1], [0 1 1 0 0]*J(i),ColorSpecs{:})
            hold on
         end
         
         % Add captions
         xlabel('x ({\mu}m)')
         ylabel('J [(A/({\mu}m)^2]')
      end
      
      function plotJfile(iew)
         % Plots exactly what in the file exported by FastHenry (j_P1.mat),
         % which may not correspond to the iew object (should maybe be a
         % static method?)
         
         filename = [iew.foldRoot,'j_P1.mat'];
         plotOn = true;
         
         J = dlmread(filename);
         J(:,1:3) = J(:,1:3)*1e6; % convert units to um
         
         if abs(std(J(:,1))/mean(J(:,1))) < 1e-6 % x is constant
            J(:,[1 4]) = [];
         elseif abs(std(J(:,2))/mean(J(:,2))) < 1e-6 % y is constant
            J(:,[2 5]) = [];
         elseif abs(std(J(:,3))/mean(J(:,3))) < 1e-6 % z is constant
            J(:,[3 6]) = [];
         end
         
         if plotOn
            switch size(J,2)
               case 4, quiver(J(:,1),J(:,2),J(:,3),J(:,4))
               case 6, quiver(J(:,1),J(:,2),J(:,4),J(:,5))
            end
            
            xlabel('x ({\mu}m)'), ylabel('y ({\mu}m)')
            hold on
         end
      end % function
      
      %% Other plotting:
      
      function plotNodes(iew)
         % Plot the locations of nodes
         plot(iew.nodes.x, iew.nodes.y,'.g')
         hold on
      end
      
      function plot(iew)
         % Plot the washer and also the current distribution
         iew.writeCIF();
         iewasher.plotCIF();
         iew.plotJ();
      end
      
      function varargout = plotContour(iew)
         % Function plots a color contour of the magnitude of the current
         % density (sqrt(Jx^2+Jy^2) by creating a 3d interpolation.
         
         % Find indices that aren't associated with crenellations:
         inds = iew.nodes.x > -iew.oDim(1)/2 & ...
            iew.nodes.y < iew.oDim(2)/2;
         % Create the interpolation:
         F = TriScatteredInterp(iew.nodes.x(inds), iew.nodes.y(inds), ...
            sqrt(iew.nodes.J(inds,1).^2 + iew.nodes.J(inds,2).^2) );
         [X,Y] = meshgrid(linspace(-iew.oDim(1),iew.oDim(1),1000),...
            linspace(-iew.oDim(2),iew.oDim(2),1000));
         % Evaluate the interpolation on the meshgrid:
         qz = F(X,Y);
         % Remove interpolated values that fall inside the hole:
         qz(abs(X) <= iew.iDim(1)/2 & abs(Y) <= iew.iDim(2)/2) = NaN;
         [~, figHan] = contourf(X, Y, qz);
         colormap cool
         axis equal
         hold on
         
         if nargout == 1, varargout{1} = figHan; end
      end
      
      function plotSegs(iew,dir,varargin)
         % Plots the segments in the discretization
      
         % Plots the horizontal or vertical segments.
         switch upper(dir)
            case 'V'
               for i = find(iew.edges.dir == 'v')'
                  plot(iew.edges.xMid(i) + ...
                     [-1 1 1 -1 -1]*iew.edges.w(i)/2, ...
                     [1 1 0 0 1]*iew.nodes.y(iew.edges.N1(i)) + ...
                     [0 0 1 1 0]*iew.nodes.y(iew.edges.N2(i)), ...
                     varargin{:})
                  hold on
               end
            case 'H'
               for i = find(iew.edges.dir == 'h')'
                  plot([1 0 0 1 1]*iew.nodes.x(iew.edges.N1(i)) + ...
                     [0 1 1 0 0]*iew.nodes.x(iew.edges.N2(i)),...
                     iew.edges.yMid(i) + ...
                     [-1 -1 1 1 -1]*iew.edges.w(i)/2, varargin{:})
                  hold on
               end
            otherwise
               disp('Direction not recognized.')
         end % switch
      end
      
      function varargout = plotStreamline(iew)
         % Function plots regularly spaced streamlines
      
         % Eliminate nodes corresponding to crenellations:
         inds = iew.nodes.x > -iew.oDim(1)/2 & iew.nodes.y < iew.oDim(2)/2;
         x = unique(iew.nodes.x(inds));
         y = unique(iew.nodes.y(inds));
         [X, Y] = meshgrid(x,y);
         Jx = nan(size(X));
         Jy = nan(size(Y));
         
         for k = find(inds)'
            i(k) = find(iew.nodes.x(k) == x);
            j(k) = find(iew.nodes.y(k) == y);
            if i(k) && j(k)
               Jx(j(k),i(k)) = iew.nodes.J(k,1);
               Jy(j(k),i(k)) = iew.nodes.J(k,2);
            end
         end
         
         for j = size(X,1):-1:1
            if all(isnan(Jx(j,:)) | Jx(j,:) == 0)
               x(j) = [];
               y(j) = [];
               X(j,:) = [];
               Y(j,:) = [];
               Jx(j,:) = [];
               Jy(j,:) = [];
            end
         end
         
         % ----- Determine starting points for streamlines -----
         Nstreams = 10;
         
         % Find indices that correspond to a slice across the bottom
         % side of the washer along x = 0 (approximately):         
         [y, J, inds] = iew.JindsAcrossLinewidth;
         s1 = interp1(y, J*iew.thickness, 'linear','pp');
         sInt = fnint(s1);
         
         % Space the streamlines by integrated current. That way, density
         % of lines corresponds to density of current.
         starty = zeros(Nstreams,1);
         for i = 1:Nstreams
            starty(i) = fzero(@(x) fnval(sInt,x) - i/(1+Nstreams), 0.5);
         end
         startx = mean(iew.edges.xMid(inds)) * ones(size(starty));
         
         % Start from the slice, then go left and right. You can do a more
         % elegant stitch using stream2, but I don't care enough to do
         % that. This works:
         streamHanL = streamline(X,Y,Jx,Jy, startx, starty);
         streamHanR = streamline(X,Y,-Jx,-Jy, startx, starty);
         streamHan = [streamHanL; streamHanR];
         set(streamHan,'color','r')
         
         if nargout == 1, varargout{1} = streamHan; end
      end
           
   end % methods
   
   methods (Static)
      %% plotCIF
      function plotCIF(filename)
         % Plot the CIF file (independent of iewasher object)
         
         if nargin == 0
            filename = 'washer.cif';
         end
         
         % Import the CIF file, ignoring line breaks and only using
         % semicolons to break statements. We'll have to remove the
         % linebreaks from the string later.
         fid = fopen(filename);
         CIFfile = textscan(fid,'%s', 'delimiter', ';', 'endofline',';');
         fclose(fid);
         CIFfile = CIFfile{1};
         
         figure
         for i = 1:length(CIFfile)
            % Remove linebreaks from the string:
            CIFfile{i} = regexprep(CIFfile{i},'\r\n','');
            
            % Search for the opening declaration of what type of object
            % the information on the line represents:
            switch upper(regexpi(CIFfile{i}, ...
                  '[^0-9]+(?= )','match','once'))
               case 'B' % box
                  % B length width x-midpoint y-midpoint;
                  vals = sscanf(CIFfile{i},'%*s %d %d %d %d');
                  if vals(2) == 0, vals(2) = 10; end
                  vals = vals/1000;
                  patch([-1 1 1 -1]*vals(1)/2 + vals(3), ...
                     [-1 -1 1 1]*vals(2)/2 + vals(4), 'c')
               case 'P' % polygon
                  % P x1 y1 x2 y2 x3 y3 ...;
                  xystr = regexp(CIFfile{i},'(-*\d+ -*\d+)|;','match');
                  xy = zeros(length(xystr),2);
                  for j = 1:length(xystr)
                     xy(j,:) = sscanf(xystr{j},'%d %d')/1000;
                  end
                  patch(xy(:,1),xy(:,2),([0,176,240]+1)/256)
            end
         end
         
         axis equal
         xlabel('x ({\mu}m)'), ylabel('y ({\mu}m)')
         hold on
         
      end
      
   end % methods (Static)
end % classdef