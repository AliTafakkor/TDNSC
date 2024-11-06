function vararginparse(args, reqargs, optargs, optdefs)
% VARARGINPARSE Parses varargin for a function, handling errors for invalid input
%
% This function processes variable input arguments (varargin) and checks 
% if they match the required and optional arguments specified. It also
% provides error handling in case of incorrect or missing inputs.
%
% INPUTS:
%   args     - Cell array of input arguments (varargin) passed to the function.
%   reqargs  - Cell array specifying the names of required arguments. The 
%              function checks if all required arguments are provided.
%   optargs  - Cell array specifying the names of optional arguments.
%   optdefs  - Cell array of default values corresponding to the optional
%              arguments. If an optional argument is not provided in args,
%              the function assigns the default value from optdefs.
%
% OUTPUTS:
%   The function does not return outputs explicitly, but it processes the
%   input arguments to make them accessible as variables within the calling 
%   function, using appropriate error handling for input validation.
%
% EXAMPLE USAGE:
%   function exampleFunction(varargin)
%       reqargs = {'input1', 'input2'};
%       optargs = {'option1', 'option2'};
%       optdefs = {10, 'defaultString'};
%       vararginparse(varargin, reqargs, optargs, optdefs);
%
%       % Now input1, input2, option1, and option2 are accessible as variables
%
% ERROR HANDLING:
%   - If any required argument is missing, an error is thrown.
%   - If an unrecognized argument is passed, an error is raised.
%   - Optional arguments not provided in args are set to their default values.
%
% Ali Tafakkor (atafakko@uwo.ca),  University of Western Ontario

% Ensure args has an even number of elements for valid name-value pairs
if(mod(length(args),2)), error('Mismatched number of arg names and values!'); end

% Extract argument names from args if not empty
if(~isempty(args)), argnames = args(1:2:end);
else, argnames = {}; end

% Parse and assign required arguments, raising an error if any are missing
if (~isempty(reqargs))
    for arg = reqargs
        ind = find(strcmp(argnames, arg{1}));
        if ind, assignin('caller', args{2*ind-1}, args{2*ind});
        else, error("Argument '%s' is required, but was not provided!", arg{1}); end
    end
end

% Parse and assign optional arguments, using default values if not provided
if (~isempty(optargs))
    for i = 1:length(optargs)
        ind = find(strcmp(argnames, optargs{i}));
        if ind, assignin('caller', args{2*ind-1}, args{2*ind});
        else, assignin('caller', optargs{i}, optdefs{i}); end
    end
end
