% function batchIEW

emailSelf = true;

try
   N = 25;
   
   %% ----- Constant W, varying R -----
%    R = [1.05 25];
%    W = 1;
%    constW1 = runBatchIEW('N',25, 'R',R,'W',W);
   
%    W = 2;
%    R = [1.05 25]*W;
%    constW2 = runBatchIEW('N',25, 'R',R,'W',W);
   
%    W = 10;
%    R = [1.05 25]*W;
%    constW10 = runBatchIEW('N',25, 'R',R,'W',W);
   
   %% ----- Constant hole, varying W -----   
%    hole = 1;
%    W = hole*[1/20, 10];
%    Ws = logspace(log10(W(1)),log10(W(2)),N);
%    Ws = round(500*Ws)/500;
%    Rs = hole/2 + Ws;
%    constHole1 = runBatchIEW('N',N, 'R',Rs, 'W',Ws);
%    
%    W = [0.1, 20];
%    Ws = logspace(log10(W(1)),log10(W(2)),N);
%    Ws = round(500*Ws)/500;
%    hole = 2;
%    Rs = hole/2 + Ws;
%    constHole2 = runBatchIEW('N',N, 'R',Rs, 'W',Ws);
%    
%    W = [0.5, 100];
%    Ws = logspace(log10(W(1)),log10(W(2)),N);
%    Ws = round(500*Ws)/500;
%    hole = 10;
%    Rs = hole/2 + Ws;
%    constHole10 = runBatchIEW('N',N, 'R',Rs, 'W',Ws);
%    
%    hole = 20;
%    W = hole*[1/20, 10];
%    Ws = logspace(log10(W(1)),log10(W(2)),N);
%    Ws = round(500*Ws)/500;
%    Rs = hole/2 + Ws;
%    constHole20 = runBatchIEW('N',N, 'R',Rs, 'W',Ws);
   
   %% ----- Constant R, varying W -----
%    R = 100;
%    W = [0.1 0.9]*R;
%    constR100 = runBatchIEW('N',25, 'R',R, 'W',W);
%    
%    R = 20;
%    W = [0.1 0.9]*R;
%    constR20 = runBatchIEW('N',25, 'R',R, 'W',W);

%    R = 10;
%    W = [0.1 0.9]*R;
%    constR10 = runBatchIEW('N',25, 'R',R, 'W',W);
%    
%    R = 2;
%    W = [0.1 0.9]*R;
%    constR2 = runBatchIEW('N',25, 'R',R, 'W',W);
   
   %% ----- Constant R/W, varying W -----
%    W = [0.5 100];
%    R = 2*W;
%    constRW2 = runBatchIEW('N',25, 'R',R, 'W',W);
%    
%    W = [0.5 100];
%    R = 20*W;
%    constRW20 = runBatchIEW('N',25, 'R',R, 'W',W);
   
   %% ----- Devices -----
   % MIT5B3:
%    W = 0.5;
%    R = [1.5 3 6 12];
%    MIT5B3 = runBatchIEW('R',R, 'W',W);
   
   %% ----- Constant R, W; varying lambda -----
   Rs = 2;
   Ws = 0.4;
   lam = [20 200]/1e3; % Don't forget that the default units are um!!!
   varyLamR2W04 = runBatchIEW('N',N, 'R',Rs, 'W',Ws, 'lam',lam); % 6ks
   
%    Rs = 3;
%    Ws = 1;
%    lam = [20 200]/1e3; % Don't forget that the default units are um!!!
%    varyLamR3W1 = runBatchIEW('N',N, 'R',Rs, 'W',Ws, 'lam',lam);
   
%    Rs = 10;
%    Ws = 3;
%    lam = [20 200]/1e3; % Don't forget that the default units are um!!!
%    varyLamR10W3 = runBatchIEW('N',N, 'R',Rs, 'W',Ws, 'lam',lam); % 8ks
   
   Rs = 45;
   Ws = 20;
   lam = [20 200]/1e3; % Don't forget that the default units are um!!!
   varyLamR45W20 = runBatchIEW('N',N, 'R',Rs, 'W',Ws, 'lam',lam); % 6ks
   
   %% ---- Constant R, W; varying b -----
   
%    Rs = 5;
%    Ws = 3;
%    b = [20 200]/1e3; % Don't forget that the default units are um!!!
%    varybR5W3 = runBatchIEW('N',N, 'R',Rs, 'W',Ws, 'b',b); % 26ks
   
   Rs = 2;
   Ws = 0.4;
   b = [20 200]/1e3; % Don't forget that the default units are um!!!
   varyb25W04 = runBatchIEW('N',N, 'R',Rs, 'W',Ws, 'b',b); % 5ks
   
catch ME
   if emailSelf
      sendGmail('santon@berkeley.edu', 'batchIEW error', ...
         {ME.identifier,ME.message});
   end
end

sendGmail('santon@berkeley.edu', 'Simulations finished!')

% end