function retval = subsref(IFSObj,S)

% Add the second argument if necessary:
if length(S) == 1
    S(2).type = '()';
    S(2).subs = {[]};
end;

% Checking validity:
if length(S) < 1 || length(S) > 2
    error('Only IndexedFaceSet_obj_name.func_name(param) is supported.');
end;

if S(1).type ~= '.' || ~ischar(S(1).subs)
    error('Only IndexedFaceSet_obj_name.func_name(param) is supported.');
end;

if ~strcmp(S(2).type,'()')
    error('Only IndexedFaceSet_obj_name.func_name(param) is supported.');
end;

try
    retval = feval(S(1).subs,IFSObj,S(2).subs{:});
catch
    error(lasterr);
end;

