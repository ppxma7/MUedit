%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decomposition of surface high-density EMG signals in motor unit discharge
% times using the method developed by Negro et al. (2016)

% FastICA: Hyvärinen & Oja (1997)
% gCKC : Holobar & Zazula (2007)

% Pseudo algorithm:
% Step 0: Load the HDsEMG data
%       0a: determine the number and type of grids
% Step 1: Preprocessing
%       1a: Signal filtering
%       1b: Removing line interference
%       1c: Differentiation in the time domain (only if the signal is good)
%       1d: Signal extension
%       1e: Removing the mean
% Step 2: Whitening
%       2a: Get eigenvalues and eigenvectors
%       2b: Zero-phase component analysis (ZCA)
% Step 3: FastICA method
%       3a: Fixed point algorithm (end when sparsness is maximized)
%       3b: Orthogonalization
%       3c: Normalization
% Step 4: Minimization of the CoV of discharge times (end when CoV is
% minimized)
%       4a: Estimate the i source
%       4b: Peak detection
%       4c: Kmean ++ classification to separate the high peaks (motor units
%       A)from small motor units (motor units B)
%       4d: Calculate CoV
%       4e: Update the separation vector
% Step 5: Select the reliable MU (SIL >0.9)
%       5a: Calculate SIL
%       5b: add the separation vector to the separation matrice
% Step 6: Update the residual signal
%       6a: Source deflation (Orthogonalization)
%       6b: go back to Step 3 to 6

%%
clear
close all;
clc;
%% Input parameters
parameters.pathname = 'C:\Users\masgh\data\'; % add a '/' at the end for Mac OS, add a '\' at the end for Windows
%parameters.filename = 'MA_190326_MMtrial_TRAP10_90DEG_2ARRAYS_3.mat'; % filename.otb+ or filename.mat
%parameters.filename = 'Injectables_002MA_VLVM_TRAP_10_TWOINJECTABLES_REALIGNED.otb4';


%multfiles = {{'rec_5_7_forDecomp.mat'; 2},{'rec_6_8_forDecomp.mat'; 2}}; % each cell contains the filename (matlab or otb) and nr of windows for this file
multfiles = {{'MA_190326_MMtrial_TRAP10_90DEG_2ARRAYS_3.mat'; 1}}; % each cell contains the filename (matlab or otb) and nr of windows for this file
nr_fil      = numel(multfiles);

run_parallel        = true; % true

if run_parallel
    desiredNumWorkers   = 2;  % or whatever number you want
end

% DECOMPOSITION PARAMETERS
parameters.NITER = 10;
parameters.ref_exist = 2; % if ref_signal exist ref_exist = 1; if not ref_exist = 0 and manual selection of windows. Michael - Add in 2, for drawrectangle version
%parameters.ref_name = 'acquired'; % MICHAEL - adding this in - actually dont need
parameters.ref_idx = 1; % 1, 4, or 7 for multichannel files (typically)
parameters.checkEMG = 0; % 0 = Consider all the channels ; 1 = Visual checking
%parameters.nwindows = 1; % number of segmented windows over each contraction
parameters.differentialmode = 0; % 0 = no; 1 = yes (filter out the smallest MU, can improve decomposition at the highest intensities
parameters.initialization = 1; % 0 = max EMG; 1 = random weights
parameters.peeloff = 0; % 0 = no; 1 = yes (update the residual EMG by removing the motor units with the highest SIL value)
parameters.covfilter = 0; % 0 = no; 1 = yes (filter out the motor units with a coefficient of variation of their ISI > than parameters.covthr)
parameters.refineMU = 0; % 0 = no; 1 = yes (refine the MU spike train over the entire signal 1-remove the discharge times that generate outliers in the discharge rate and 2- reevaluate the MU pulse train)
parameters.drawingmode = 0; % 0 = Output in the command window ; 1 = Output in a figure
parameters.duplicatesbgrids = 0; % 0 = do not consider duplicates between grids ; 1 = Remove duplicates between grids
parameters.alignMUAP = 0; % 0 = no; 1 = yes (realign the discharge time with the peak of the MUAP (channel with the MUAP with the highest p2p amplitude from double diff EMG signal)

% SPECIFIC VALUES
parameters.thresholdtarget = 0.8; % threshold to segment the target displayed to the participant, 1 being the maxima of the target (e.g., plateau)
parameters.nbextchan = 1000; % nb of extended channels (1000 in Negro 2016, can be higher to improve the decomposition)
parameters.edges = 0.2; % edges of the signal to remove after preprocessing the signal (in sec)
parameters.contrastfunc = 'skew'; % contrast functions: 'skew', 'kurtosis', 'logcosh'
parameters.silthr = 0.90; % Threshold for SIL values
parameters.covthr = 0.5; % Threshold for CoV of ISI values
parameters.peeloffwin = 0.025; % duration of the window (ms) for detecting the action potentials from the EMG signal
parameters.duplicatesthresh = 0.3; % threshold that define the minimal percentage of common discharge times between duplicated motor units
parameters.CoVDR = 0.3; % threshold that define the CoV of Discharge rate that we want to reach for cleaning the MU discharge times when refineMU is on


alldata = cell(nr_fil,1);
for fil = 1:nr_fil

    parameters.filename = multfiles{fil}{1};
    parameters.nwindows = multfiles{fil}{2};

    f = waitbar(0,'Preprocessing - Load the HDsEMG data');


    %% Step 0: Load the HDsEMG data
    %       0a: determine the number and type of grids
    C = strsplit(parameters.filename,'.');
    if isequal(C{end}, 'mat')
        load([parameters.pathname parameters.filename], 'signal');
    elseif isequal(C{end}, 'otb+')
        [~, signal] = openOTBplus(parameters.pathname, parameters.filename,0);
    elseif isequal(C{end},'otb4')
        [~, signal] = OpenOTB4_nottingham(parameters.pathname, parameters.filename, 0);
    else
        error('cant parse the file extension!')
    end

    [signal.coordinates, signal.IED, signal.EMGmask, signal.emgtype] = formatsignalHDEMG(signal.data, signal.gridname, signal.fsamp, parameters.checkEMG);
    

    parameters.nbelectrodes = signal.nChan / signal.ngrid;
    channelstore = reshape(1:parameters.nbelectrodes*signal.ngrid, [], signal.ngrid)';
    parameters.channels = mat2cell(channelstore, ones(1, signal.ngrid));

    arraynb = zeros(size(signal.data,1),1);
    ch1 = 1;
    for i = 1:signal.ngrid
        arraynb(ch1:ch1+length(signal.EMGmask{i})-1) = i;
        ch1 = ch1+length(signal.EMGmask{i});
    end
    %%
    % Step 0 <opt> Selection of the region of interest
    if parameters.ref_exist == 1
        signalprocess.ref_signal = signal.target;
        signalprocess.coordinatesplateau = segmenttargets(signalprocess.ref_signal, parameters.nwindows, parameters.thresholdtarget);
        for nwin = 1:length(signalprocess.coordinatesplateau)/2
            for i = 1:signal.ngrid
                signalprocess.data{i,nwin} = signal.data(arraynb==i, signalprocess.coordinatesplateau(nwin*2-1):signalprocess.coordinatesplateau(nwin*2));
                signalprocess.data{i,nwin}(signal.EMGmask{i} == 1,:) = [];
            end
        end
    elseif parameters.ref_exist == 2
        %idx = find(contains(signal.auxiliaryname, parameters.ref_name), 1);
        % be defensive here in case user picks value greater than number of
        % targets
        assert(parameters.ref_idx <= length(signal.auxiliaryname), ...
            'ref_idx (%d) exceeds number of auxiliary signals (%d)', ...
            parameters.ref_idx, length(signal.auxiliaryname));
        idx = parameters.ref_idx;
        signalprocess.ref_signal = signal.auxiliary(idx,:);
        signalprocess.ref_signal = (signalprocess.ref_signal - signalprocess.ref_signal(1)) / ...
            max(signalprocess.ref_signal - signalprocess.ref_signal(1));
        signal.target = signalprocess.ref_signal;
        signal.path = signalprocess.ref_signal;

        figure;
        plot(signalprocess.ref_signal, 'Color', [0.5 0.5 0.5], 'LineWidth', 2)
        ylim([0 max(signalprocess.ref_signal)*1.2])
        title([signal.auxiliaryname{idx} ' - Select windows'])
        grid on

        signalprocess.coordinatesplateau = zeros(1, parameters.nwindows*2);
        for nwin = 1:parameters.nwindows
            title([signal.auxiliaryname{idx} ' - Select window #' num2str(nwin)])
            roi = drawrectangle('DrawingArea', 'auto');
            x = sort([roi.Position(1), roi.Position(1) + roi.Position(3)], 'ascend');
            x(x<1) = 1;
            x(x>length(signalprocess.ref_signal)) = length(signalprocess.ref_signal);
            signalprocess.coordinatesplateau(nwin*2-1) = floor(x(1));
            signalprocess.coordinatesplateau(nwin*2) = floor(x(2));
            for i = 1:signal.ngrid
                signalprocess.data{i,nwin} = signal.data(arraynb==i, signalprocess.coordinatesplateau(nwin*2-1):signalprocess.coordinatesplateau(nwin*2));
                signalprocess.data{i,nwin}(signal.EMGmask{i} == 1,:) = [];
            end
        end
        clearvars x idx
    else
        tmp = zeros(floor(size(signal.data,1)/2), size(signal.data,2));
        for i = 1:floor(size(signal.data,1)/2)
            tmp(i,:) = movmean(abs(signal.data(i,:)), signal.fsamp);
        end
        signalprocess.ref_signal = mean(tmp,1);
        signal.target = signalprocess.ref_signal;
        signal.path = signalprocess.ref_signal;
        plot(tmp','Color', [0.5 0.5 0.5],'LineWidth',0.5)
        hold on
        plot(signalprocess.ref_signal,'k','LineWidth',1)
        ylim([0 max(signalprocess.ref_signal)*1.5])
        title('EMG amplitude for 50% of the EMG channels')
        grid on
        clearvars tmp maskEMG
        [A,B] = ginput(2*parameters.nwindows);
        A(A<1) = 1;
        A(A>length(signalprocess.ref_signal)) = length(signalprocess.ref_signal);
        signalprocess.coordinatesplateau = zeros(1,parameters.nwindows*2);
        for nwin = 1:parameters.nwindows
            signalprocess.coordinatesplateau(nwin*2-1) = floor(A(nwin*2-1));
            signalprocess.coordinatesplateau(nwin*2) = floor(A(nwin*2));
            for i = 1:signal.ngrid
                signalprocess.data{i,nwin} = signal.data(arraynb==i, signalprocess.coordinatesplateau(nwin*2-1):signalprocess.coordinatesplateau(nwin*2));
                signalprocess.data{i,nwin}(signal.EMGmask{i} == 1,:) = [];
            end
        end
        clearvars A B
    end
    
    alldata{fil}.signal = signal;
    alldata{fil}.signalprocess = signalprocess;
    clear signal signalprocess;
    close all hidden;

end

clear i nwin f
keepvariables = who;
keepvariables{end+1} = 'keepvariables';


%%

for fil = 1:nr_fil
    clearvars('-except', keepvariables{:});
    parameters.filename = multfiles{fil}{1};
    parameters.nwindows = multfiles{fil}{2};
    signal              = alldata{fil}.signal;
    signalprocess       = alldata{fil}.signalprocess;
    alldata{fil} = [];
    disp(parameters.filename);
    tic
        PulseT_grd      = cell(1,signal.ngrid);
        distimenew_grd  = cell(signal.ngrid,1);
        if run_parallel
            % Check if a pool is already active
            p = gcp('nocreate');  % Don't create a pool if none exists
            
            if isempty(p)
                % No pool, create one with desired number of workers
                parpool('local', desiredNumWorkers);
            elseif p.NumWorkers ~= desiredNumWorkers
                % A pool exists, but with the wrong number of workers — optional:
                delete(p);
                parpool('local', desiredNumWorkers);
            end
        
            % Shared constants (only needed once per worker)
            signalConst = parallel.pool.Constant(@() signal);
            paramConst  = parallel.pool.Constant(@() parameters);
            
            % Prepare cell array of per-grid signalprocess structs
            signalprocesscell = cell(1, signal.ngrid);
            for grid_nr = 1:signal.ngrid
                sp = signalprocess;
                for grd = 1:signal.ngrid
                    if grd ~= grid_nr
                        sp.data{grd} = [];  % Clear all other grid entries
                    end
                end
                signalprocesscell{grid_nr} = sp;
            end
            clear sp grd grid_nr

            parfor grid_nr = 1:signal.ngrid
                [PulseT_grd{grid_nr}, distimenew_grd{grid_nr}] = decompGrid(signalConst.Value, paramConst.Value, signalprocesscell{grid_nr}, grid_nr, run_parallel);
            end % end of "for each grid"
            clear signalConst paramConst signalprocesscell
        else
            for grid_nr = 1:signal.ngrid
                [PulseT_grd{grid_nr}, distimenew_grd{grid_nr}] = decompGrid(signal, parameters, signalprocess, grid_nr, run_parallel);
            end % end of "for each grid"
        end%if-->run_parallel
        toc
        signal.Pulsetrain = PulseT_grd;
        maxnumber = max(cellfun(@(dist_grd) numel(dist_grd), distimenew_grd));
        for grd_nr = 1:signal.ngrid
            signal.Dischargetimes(grd_nr,:) = [distimenew_grd{grd_nr}, cell(1, maxnumber-numel(distimenew_grd{grd_nr}))];
        end

        savename = fullfile(parameters.pathname, [parameters.filename '_decomp.mat']);
        save(savename, 'signal', 'parameters', '-v7.3');

end

