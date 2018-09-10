function [ch] = getkeySH
%aka getkeyrozzo

% check the input arguments
callstr = 'set(gcbf,''Userdata'',double(get(gcbf,''Currentcharacter''))) ; uiresume ' ;

% Set up the figure
fh = figure(...
    'name','Press a key', ...
    'keypressfcn',callstr, ...
    'windowstyle','modal',...
    'numbertitle','off', ...
    'position',[0 0  1 1],...
    'userdata','timeout') ;

% Wait for a key press before uiresume is  executed
uiwait;
ch = get(fh,'Userdata') ;  % and the key itself

% clean up the figure, if it still exists
if ishandle(fh);    delete(fh); end
end