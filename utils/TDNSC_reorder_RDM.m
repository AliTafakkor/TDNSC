function [RDM_reord, labels_reord] = TDNSC_reorder_RDM(p, RDM, labels, indmap, level)
% TDNSC_REORDER_RDM changes the order of conditions / stimuli in a RDM. 
%
% INPUTS:
%   p      (required) - Struct specifying the project parameters.
%   RDM    (required) - Square matrix representating the RDM.
%   labels (required) - Vector or cell array specifying the labels of RDM
%                       rows/columns.
%
%   indmap (optional) - Vector specifying mapping of the indecies 
%                       (must be the same size as the RDM rows/columns).
%                       see NOTES for the default behavior.
%   level  (optional) - String specifying the decoding level 
%                       function checks if all required arguments
%                       are provided.
%
% OUTPUTS:
%   RDM_reord         - RDM with the new order of condition / stimuli.
%   labels_reord      - New order of condition / stimuli labels.
%
% EXAMPLE USAGE:
%   
%
% NOTE:
%   - If indmap and level are not given, level will be automatically set
%     according to the RDM size.
%   - Default map will move the last 20 stimuli / category after the first
%     20 stimuli / category to change the category order from "Animals,
%     Objects, Scenes, People" to "Animals, People, Objects, Scenes". The
%     first order was used in conducting the experience, but later it was
%     decided that it's better to have "animate" categories next to each
%     other.
%
% Ali Tafakkor (atafakko@uwo.ca),  University of Western Ontario

% Check inputs
if nargin < 2, error("Arguments 'p' and 'RDM' are required!"); end

n = size(RDM,1);
if size(RDM,2) ~= n, error('RDM must be an square matrix!'); end

if nargin < 5
    % Figure out level automatically
    if n == p.stim.num_cat, level='category';
    elseif n == p.stim.num, level='stimuli';
    else, error("Decoding level cannot be found!"); end
end

if nargin < 4
    % Set default indmap default to reorder
    % map to move people before objects and scenes,
    % new order (Animals,People,Objects,Scenes)
    if strcmp(level, 'category')
        indmap = [1 4 2 3];
    elseif strcmp(level, 'stimuli')
    % map to move people before objects and scenes
    % new order (Animals,People,Objects,Scenes)
        indmap = [1:20 61:80 21:60];
    end
end

% Apply map
RDM_reord = RDM(indmap, indmap);

if nargin > 2, labels_reord = labels(indmap); end