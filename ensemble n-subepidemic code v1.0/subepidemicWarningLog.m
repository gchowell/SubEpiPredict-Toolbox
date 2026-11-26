function varargout = subepidemicWarningLog(action,varargin)
%SUBEPIDEMICWARNINGLOG Keep the toolbox's warnings of one run and print them at its end.
%   subepidemicWarningLog('start')            empties the list (start of a driver or script)
%   subepidemicWarningLog('add',id,message)   called by subepidemicWarning.m
%   subepidemicWarningLog('summary',label)    prints the warnings raised since 'start', one
%                                             line per kind with the number of times and
%                                             the first message; prints nothing if none
%   L = subepidemicWarningLog('get')          returns the list (struct array: id, message, count)
%
%   A run prints a lot between its first line and its last, and a warning
%   raised in the middle is easy to miss. The fit driver and the four scripts
%   call 'start' at their beginning and 'summary' at their end. Only warnings
%   raised through subepidemicWarning.m are listed (those with an identifier
%   starting with SubEpiPredict:), not MATLAB's own.

persistent entries
if isempty(entries)
    entries = struct('id',{},'message',{},'count',{});
end
switch action
    case 'start'
        entries = struct('id',{},'message',{},'count',{});
    case 'add'
        id = varargin{1};
        k = find(strcmp({entries.id},id),1);
        if isempty(k)
            entries(end+1) = struct('id',id,'message',varargin{2},'count',1);
        else
            entries(k).count = entries(k).count+1;
        end
    case 'summary'
        if nargin>1, label = varargin{1}; else, label = 'SubEpiPredict'; end
        if ~isempty(entries)
            fprintf('\n%s finished with %d kind(s) of warning:\n',label,numel(entries));
            for k = 1:numel(entries)
                if entries(k).count==1
                    fprintf('  %s: %s\n',entries(k).id,entries(k).message);
                else
                    fprintf('  %s (%d times; first:) %s\n',entries(k).id,entries(k).count,entries(k).message);
                end
            end
        end
    case 'get'
        varargout{1} = entries;
    otherwise
        error('SubEpiPredict:WarningLog','Unknown action "%s".',action);
end
end
