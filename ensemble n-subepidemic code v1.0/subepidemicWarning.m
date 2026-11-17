function subepidemicWarning(id,messageFormat,varargin)
%SUBEPIDEMICWARNING Raise a warning of the toolbox and remember it for the summary.
%   subepidemicWarning(id,messageFormat,A1,...) is warning(id,messageFormat,A1,...)
%   and also records the warning with subepidemicWarningLog.m, so that the
%   driver or script can list it again at its end. A warning switched off
%   with warning('off',id) is neither shown nor recorded.

message = sprintf(messageFormat,varargin{:});
state = warning('query',id);
if strcmp(state.state,'on')
    subepidemicWarningLog('add',id,message);
end
warning(id,'%s',message);
end
