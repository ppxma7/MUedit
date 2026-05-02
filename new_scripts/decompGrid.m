function [PulseT_grd, distimenew_grd] = decompGrid(signal, parameters, signalprocess, grid_nr, parallel)

disp(['grid: ', num2str(grid_nr)]);
for nwin = 1:length(signalprocess.coordinatesplateau)/2

    % Step 1: Preprocessing
    %       1a: Removing line interference (Notch filter)
    if ~parallel
        f = waitbar(0.2, ['Grid #' num2str(grid_nr) ' - Preprocessing - Filtering the HDsEMG data']);
    end%if parallel
        signalprocess.data{grid_nr,nwin} = notchsignals(signalprocess.data{grid_nr,nwin},signal.fsamp);
    %       1b: Bandpass filtering
        signalprocess.data{grid_nr,nwin} = bandpassingals(signalprocess.data{grid_nr,nwin},signal.fsamp, signal.emgtype(grid_nr));
    
    %       1c: Differentiation (perform only if there is many motor units,
    %      filter out the smallest motor units) useful for high intensities
        if parameters.differentialmode == 1
                signalprocess.data{grid_nr,nwin} = diff(signalprocess.data{grid_nr,nwin},1,2);
        end
    
    %       1d: Signal extension (extension factor calculated to reach 1000
    %       channels)
    if ~parallel
        waitbar(0.4, f, ['Grid #' num2str(grid_nr) ' - Preprocessing - Extending the HDsEMG data'])
    end%if-->parallel
    
        signalprocess.exFactor = round(parameters.nbextchan/size(signalprocess.data{grid_nr,nwin},1));
        signalprocess.ReSIG = zeros(signalprocess.exFactor * size(signalprocess.data{grid_nr,nwin},1));
        signalprocess.iReSIG{nwin} = zeros(signalprocess.exFactor * size(signalprocess.data{grid_nr,nwin},1));
        
        signalprocess.eSIG{nwin} = extend(signalprocess.data{grid_nr,nwin},signalprocess.exFactor);
        signalprocess.ReSIG = signalprocess.eSIG{nwin} * signalprocess.eSIG{nwin}' / size(signalprocess.eSIG{nwin},2);
        signalprocess.iReSIG{nwin} = pinv(signalprocess.ReSIG);
    
    %       1e: Removing the mean
        signalprocess.eSIG{nwin} = demean(signalprocess.eSIG{nwin});
    
    % Step 2: Whitening
    %Whitening with a regularization factor (average of the smallest half of 
    %the eigenvalues of the covariance matrix from the extended signals)
    
    %       2a: Get eigenvalues and eigenvectors (regularization factor =>
    %       average smallest half of eigenvalues)
    if ~parallel
        waitbar(0.6, f, ['Grid #' num2str(grid_nr) ' - Preprocessing - Whitening the HDsEMG data'])
    end%if-->parallel
    [E, D] = pcaesig(signalprocess.eSIG{nwin}); %Returns the eigenvector (E) and diagonal eigenvalue (D) matrices
    
    %       2b: Zero-phase component analysis
        [signalprocess.wSIG{nwin}, ~, ~] = whiteesig(signalprocess.eSIG{nwin}, E, D);
        clearvars E D
    
    % Remove the edges
        signalprocess.eSIG{nwin} = signalprocess.eSIG{nwin}(:,round(signal.fsamp*parameters.edges):end-round(signal.fsamp*parameters.edges));
        signalprocess.wSIG{nwin} = signalprocess.wSIG{nwin}(:,round(signal.fsamp*parameters.edges):end-round(signal.fsamp*parameters.edges));
    
    % if grid_nr == 1
        signalprocess.coordinatesplateau(nwin*2-1) = signalprocess.coordinatesplateau(nwin*2-1) + round(signal.fsamp*parameters.edges)-1;
        signalprocess.coordinatesplateau(nwin*2) = signalprocess.coordinatesplateau(nwin*2) - round(signal.fsamp*parameters.edges);
    % end
   
    
    % Step 3: FastICA method
    if parallel
        disp(['Grid ', num2str(grid_nr), ' window ', num2str(nwin), ' starting decomposition'])
    else
        waitbar(0.8, f, ['Grid #' num2str(grid_nr) ' - Decomposition - Decomposing the HDsEMG data'])
    end%if
    
    % Initialize matrix B (n x m) n: separation vectors, m: iterations 
    % Initialize matrix MUFilters to only save the reliable filters
    % Intialize SIL and PNR
    
    signalprocess.B = zeros(size(signalprocess.wSIG{nwin},1), parameters.NITER); % all separation vectors
    signalprocess.MUFilters{nwin} = zeros(size(signalprocess.wSIG{nwin},1), parameters.NITER); % only reliable vectors
    signalprocess.w = zeros(size(signalprocess.wSIG{nwin},1), 1);
    signalprocess.icasig = zeros(parameters.NITER, size(signalprocess.wSIG{nwin},2));
    signalprocess.SIL{nwin} = zeros(1, parameters.NITER);
    signalprocess.CoV{nwin} = zeros(1, parameters.NITER);
    idx1 = zeros(1, parameters.NITER);
    
    % Find the index where the square of the summed whitened vectors is
    % maximized and initialize W with the whitened observations at this time
    if exist('f','var')
        close (f)
    end%if
    
    if parameters.drawingmode == 0 && ~parallel
        f = waitbar(0,['Grid #' num2str(grid_nr) ' - Decomposition']);
        xwb = 1/parameters.NITER;
    end
    
    
    for j = 1:parameters.NITER
        if j == 1
            signalprocess.X = signalprocess.wSIG{nwin}; % Initialize X (whitened signal), then X: residual
            if parameters.initialization == 0
                actind = sum(signalprocess.X,1).^2;
                [~, idx1(j)] = max(actind);
                signalprocess.w = signalprocess.X(:, idx1(j)); % Initialize w
            else
                signalprocess.w = randn(size(signalprocess.X,1),1); % Initialize w
            end
            time = linspace(0,size(signalprocess.X,2)/signal.fsamp,size(signalprocess.X,2));
        else
            if parameters.initialization == 0
                actind(idx1(j-1)) = 0; % remove the previous vector
                [~, idx1(j)] = max(actind);
                signalprocess.w = signalprocess.X(:, idx1(j)); % Initialize w
            else
                signalprocess.w = randn(size(signalprocess.X,1),1); % Initialize w
            end
        end
        
        signalprocess.w = signalprocess.w - signalprocess.B * signalprocess.B' * signalprocess.w; % Orthogonalization
        signalprocess.w = signalprocess.w / norm(signalprocess.w); % Normalization
        
        %       3a: Fixed point algorithm (end when sparsness is maximized)
        maxiter = 500; % max number of iterations for the fixed point algorithm
        signalprocess.w = fixedpointalg(signalprocess.w, signalprocess.X, signalprocess.B , maxiter, parameters.contrastfunc);
        
        % Step 4: Minimization of the CoV of discharge times (end when CoV is
        % minimized)
        
        % Initialize CoV (variation of interspike intervals, %) Step 4a => 4e
        [signalprocess.icasig, signalprocess.spikes] = getspikes(signalprocess.w, signalprocess.X, signal.fsamp);
        
        if length(signalprocess.spikes) > 10
            ISI = diff(signalprocess.spikes/signal.fsamp); % Interspike interval
            signalprocess.CoV{nwin}(j) = std(ISI)/mean(ISI); % Coefficient of variation
            Wini = sum(signalprocess.X(:,signalprocess.spikes),2); % update W by summing the spikes
                           
                    % Minimization of the CoV of discharge times (end when CoV is
            % minimized)
            [signalprocess.MUFilters{nwin}(:,j), signalprocess.spikes, signalprocess.CoV{nwin}(j)] = minimizeCOVISI(Wini, signalprocess.X, signalprocess.CoV{nwin}(j), signal.fsamp);
            signalprocess.B(:,j) = signalprocess.w;
        
            % Calculate SIL values
            [signalprocess.icasig, signalprocess.spikes, signalprocess.SIL{nwin}(j)] = calcSIL(signalprocess.X, signalprocess.MUFilters{nwin}(:,j), signal.fsamp);
            
            % Peel-off of the (reliable) source    
            if parameters.peeloff == 1 && signalprocess.SIL{nwin}(j) > parameters.silthr
               signalprocess.X = peeloff(signalprocess.X, signalprocess.spikes, signal.fsamp, parameters.peeloffwin);
            end
        
            if parameters.drawingmode == 1 && ~parallel
                subplot(2,1,1)
                plot(signal.target, 'k--', 'LineWidth', 2)
                line([signalprocess.coordinatesplateau(nwin*2-1) signalprocess.coordinatesplateau(nwin*2-1)],[0 max(signal.target)], 'Color', 'r', 'LineWidth', 2)
                line([signalprocess.coordinatesplateau(nwin*2) signalprocess.coordinatesplateau(nwin*2)],[0 max(signal.target)], 'Color', 'r', 'LineWidth', 2)
                title(['Grid #' num2str(grid_nr) ' - Iteration #' num2str(j) ' - Sil = ' num2str(signalprocess.SIL{nwin}(j)) ' CoV = ' num2str(signalprocess.CoV{nwin}(j))]);
                subplot(2,1,2)
                plot(time,signalprocess.icasig,time(signalprocess.spikes),signalprocess.icasig(signalprocess.spikes),'o');
                drawnow;
            elseif ~parallel
                waitbar(xwb*j, f, ['Grid #' num2str(grid_nr) ' - Iteration #' num2str(j) ' - Sil = ' num2str(signalprocess.SIL{nwin}(j)) ' CoV = ' num2str(signalprocess.CoV{nwin}(j))])
            end
        else
            signalprocess.B(:,j) = signalprocess.w;
        end
        if parallel && any([round(parameters.NITER/4), round(parameters.NITER/2), round(parameters.NITER*3/4)] == j)
            disp(['Grid ' num2str(grid_nr) ' - Iteration ' num2str(j)])
        end%if

    end%for-->j
    
    if parameters.drawingmode == 0 && exist('f', 'var')
        close(f);
    end
    
    
    
    % Filter out MUfilters below the SIL threshold
    signalprocess.MUFilters{nwin}(:,signalprocess.SIL{nwin} < parameters.silthr) = [];
    % disp(find(signalprocess.SIL{nwin} >= parameters.silthr));
    if parameters.covfilter == 1
        signalprocess.CoV{nwin}(signalprocess.SIL{nwin} < parameters.silthr) = [];
        signalprocess.MUFilters{nwin}(:,signalprocess.CoV{nwin} > parameters.covthr) = [];
    end

    if length(signalprocess.coordinatesplateau)/2>1
        thisPlateau = [signalprocess.coordinatesplateau(nwin*2-1),signalprocess.coordinatesplateau(nwin*2)];
        [PulseT, distime] = batchprocessfilters_penalty(signalprocess.MUFilters(nwin), signalprocess.wSIG(nwin), thisPlateau, signalprocess.exFactor, parameters.differentialmode, size(signal.data,2), signal.fsamp);
        
        if ~isempty(PulseT)
            [PulseT, distimenew, survivor_ind, unique_ind] = remduplicates_track(PulseT, distime, distime, round(signal.fsamp/40), 0.00025, parameters.duplicatesthresh, signal.fsamp, 'cov');
        else
            PulseT      = [];
            distimenew  = [];
        end
        signalprocess.MUFilters{nwin}=signalprocess.MUFilters{nwin}(:,[survivor_ind,unique_ind]);
        sil = zeros(size(PulseT,1),1);
        for t = 1:size(PulseT,1)
            sil(t) = getsil_MUedit_sumdpenalty(PulseT(t,:), signal.fsamp);
        end%for-->t
        goodsil = find(sil>=parameters.silthr);
        PulseT  = PulseT(goodsil,:);
        distimenew = distimenew(goodsil);
        signalprocess.MUFilters{nwin}=signalprocess.MUFilters{nwin}(:,goodsil);
        % distimenew = {distimenew};
        % distimenew = cellfun(@(dtn) dtn-thisPlateau(1), distimenew, 'Uniformoutput', false);
        % signalprocess.MUFilters{nwin} = getMUfilters_MUedit(signalprocess.data{grid_nr,nwin}, distimenew, parameters.nbextchan);
        % for dit = 1:numel(distimenew)
        %     [~, ~, silnew(dit)] = calcSIL(signalprocess.X, signalprocess.MUFilters{nwin}(:,dit), signal.fsamp);
        % end%for..>dit
        clear PulseT distimenew 
    end%if




end%for-->nwin
if ~parallel
    f = waitbar(0.8,['Grid #' num2str(grid_nr) ' - Postprocessing']);
else
    disp(['Grid #' num2str(grid_nr) ' - Postprocessing'])
end%if

% Batch processing over each window
% [PulseT, distime] = batchprocessfilters(signalprocess.dewhiteningMatrix, signalprocess.MUFilters, signalprocess.iReSIG, signalprocess.eSIG, signalprocess.coordinatesplateau, signalprocess.exFactor, parameters.differentialmode, size(signal.data,2), signal.fsamp);
[PulseT, distime] = batchprocessfilters_penalty(signalprocess.MUFilters, signalprocess.wSIG, signalprocess.coordinatesplateau, signalprocess.exFactor, parameters.differentialmode, size(signal.data,2), signal.fsamp);
% after Batch processing calc sil again and exclude to low ones --> paul do this if there where more than 1 windows
% if length(signalprocess.coordinatesplateau)/2>1
%     sil = zeros(size(PulseT,1),1);
    % for t = 1:size(PulseT,1)
    %     sil(t) = getsil_MUedit(nonzeros(PulseT(t,:))', signal.fsamp);
    % end%for-->t
%     goodsil = find(sil>=parameters.silthr);
%     PulseT  = PulseT(goodsil,:);
%     distime = distime(goodsil);
% end%if

% Realign the discharge times with the center of the MUAP
if parameters.alignMUAP == 1
    % distimea = realignspikes(signal.data((grid_nr-1)*parameters.nbelectrodes+1:(grid_nr-1)*parameters.nbelectrodes+length(signal.EMGmask{grid_nr}), :), signal.EMGmask{grid_nr}, signal.coordinates{grid_nr}, distime, signal.fsamp, parameters.peeloffwin);
    distimea = realignspikes(signal.data(parameters.channels{grid_nr}, :), signal.EMGmask{grid_nr}, signal.coordinates{grid_nr}, distime, signal.fsamp, parameters.peeloffwin); % paul changed
else
    distimea = distime;
end

% Remove duplicates
if ~isempty(PulseT)
    [PulseT, distimenew] = remduplicates(PulseT, distime, distimea, round(signal.fsamp/40), 0.00025, parameters.duplicatesthresh, signal.fsamp);
else
    PulseT      = [];
    distimenew  = [];
end

if parameters.refineMU == 1    
    % Remove outliers generating irrelevant discharge rates before manual
    % edition (1st time)
    distimenew = remoutliers(PulseT, distimenew, parameters.CoVDR, signal.fsamp);

    % Reevaluate all the unique motor units over the contractions
    % [signal.Pulsetrain{grid_nr}, distimenew] = refineMUs(signal.data((grid_nr-1)*parameters.nbelectrodes+1:(grid_nr-1)*parameters.nbelectrodes+length(signal.EMGmask{grid_nr}), :), signal.EMGmask{grid_nr}, PulseT, distimenew, signal.fsamp);
    [s.Pulsetrain{grid_nr}, distimenew] = refineMUs(signal.data(parameters.channels{grid_nr}, :), signal.EMGmask{grid_nr}, PulseT, distimenew, signal.fsamp); % paul changed
    
    % Remove outliers generating irrelevant discharge rates before manual
    % edition (2nd time)
    distimenew = remoutliers(s.Pulsetrain{grid_nr}, distimenew, parameters.CoVDR, signal.fsamp);
else
    s.Pulsetrain{grid_nr} = PulseT;
end

% Save the results
if exist("distimenew","var")
    for j = 1:length(distimenew)
        s.Dischargetimes{grid_nr,j} = distimenew{j};
    end
end
if isfield(s, 'Dischargetimes')
    distimenew_grd  = s.Dischargetimes(grid_nr,:);
else
    distimenew_grd  = [];
end%if

PulseT_grd      = s.Pulsetrain{grid_nr};

end%function