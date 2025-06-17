function [dToff, dTon, tStart, tEnd] = makeTs(t, linkIdx)
% calculates time-related parameters of binding events in a cluster based
% on the number of blink events and their link indices  

% INPUT
% t: frame number of blink
% linkIdx: indices of blinks that belong to the same binding event

% OUTPUT
% dToff: Times between binding events
% dTon: Duration of a binding event
% tStart: Start frame of a binding event
% tEnd: End frame of a binding event

% get all the unique blink IDs
blinkIDs = unique(linkIdx);

% get number of unique blink IDs
nIDs = numel(blinkIDs);

% setup arrays to hold the start and end time for each binding event
tStart = zeros(nIDs, 1);
tEnd = zeros(nIDs, 1);

% loop over all unique blink IDs
for idt = 1:nIDs

    tmpIdx = linkIdx == blinkIDs(idt);
    blinkT = t(tmpIdx);
    tStart(idt) = blinkT(1);
    tEnd(idt) = blinkT(end);

end

tSeq = [tStart'; tEnd'];

dts = diff(tSeq(:));
dTon = dts(1:2:end) + 1;

% find the off time between binding events
dToff = dts(2:2:end) - 1;
