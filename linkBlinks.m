function linkIdx = linkBlinks(t, nGap)
% This function groups blinks in a cluster into binding events based on a
% threshold "nGap" of dark frames between blink events

% INPUT
% t: frame number
% nGap: number of dark frames between blinks

% OUTPUT:
% linkIdx: indices of blinks for each binding event

% ensure that a cluster doesn't have two blinks happening at the same time
% (temporal overlap)
if isequal(size(t), size(unique(t)))

    % initialize a column with zeros for the link index
    linkIdx = zeros(size(t));

    % initialize the link group tag
    groupID = 0;

    % loop over all blinks until they're all assigned to a link index
    while ~all(linkIdx)

        % find the index of the first blink that does not have a link index
        aIdx = find(~linkIdx, 1, 'first');

        % asign a new ID
        groupID = groupID + 1;

        % set the link Index to the new ID
        linkIdx(aIdx)=groupID;

        % start the gap counter
        gapCounter=0;

        % set the current Frame to the selected blink
        cFC = t(aIdx);

        % stop when the gap becomes bigger than nGap
        while gapCounter <= nGap

            % go to the next frame
            cFC = cFC + 1;

            % find the blinks on that frame
            bIdx = find(t == cFC);

            if isempty(bIdx)
                % if there are not any blinks increase the gap counter we have a gap frame
                gapCounter=gapCounter+1;
            else
                % if there is blink reset the gap counter to 0
                gapCounter = 0;

                % assign the same group ID to the link Index of this blink
                linkIdx(bIdx) = groupID;
            end

        end

    end
else

    % case for temporal overalp (two binks happen at the same time in one
    % cluster
    linkIdx = [];

    if isempty(t)
        fprintf('qPAINTlinkBlinks failed empty T array.\n')
    else
        fprintf('qPAINTlinkBlinks failed due to temporal overlap of blinks.\n')
    end

end