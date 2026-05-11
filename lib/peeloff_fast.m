%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peel off of the Motor unit spike train

% Input: 
%   X = whitened signal
%   spikes = discharge times of the motor unit
%   fsamp = sampling frequency
%   win = window to identify the motor unit action potential with spike trigger averaging


% Output:
%   X = residual of the whitened signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = peeloff_fast(X, spikes, fsamp, win)

% parameters
    windowl = round(win * fsamp);
    L = 2 * windowl + 1;
    [nCh, nSamples] = size(X);

    % remove invalid spikes
    spikes = spikes(spikes > windowl & spikes <= nSamples - windowl);
    if isempty(spikes)
        return
    end

    % ---- Build index matrix (no loops) ----
    idx = spikes(:) + (-windowl:windowl);   % [nSpikes x L]

    % ---- Extract MUAPs for ALL channels at once ----
    % Result: [nCh x nSpikes x L]
    MUAPs = reshape(X(:, idx(:)), nCh, length(spikes), L);

    % ---- Spike-triggered average ----
    waveform = squeeze(mean(MUAPs, 2));     % [nCh x L]
    
    % ---- Reconstruct and peel off ----
    EMGtemp = zeros(size(X));
    % for ch = 1:nCh
    %     EMGtemp(ch,:) = conv(firings, waveform(ch,:), 'same');
    % end
    for k = 1:size(idx,1)
        EMGtemp(:, idx(k,:)) = EMGtemp(:, idx(k,:)) + waveform;
    end
    
    X = X - EMGtemp;
end%function
