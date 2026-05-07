function pass = test_redbird(testname, fhandle, expected, varargin)
%
% pass = test_redbird(testname, fhandle, expected, arg1, arg2, ...)
%
% Tiny test helper used by run_redbird_test.
%
% Calls fhandle(arg1, arg2, ...) and compares its first return value to
% `expected`. Numerical comparisons use a relative + absolute tolerance.
%
% Special expected values:
%     'noerror' : just ensure the call returns without throwing
%     'error'   : ensure the call DOES throw
%
% Tolerances default to 1e-9 (rel) and 1e-12 (abs); override by setting
% globals RB_RTOL / RB_ATOL before invoking.
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% license:
%     GPL version 3, see LICENSE_GPLv3.txt files for details
%
% -- this function is part of Redbird-m toolbox
%

global RB_RTOL RB_ATOL RB_FAIL RB_TOTAL
if isempty(RB_RTOL)
    RB_RTOL = 1e-9;
end
if isempty(RB_ATOL)
    RB_ATOL = 1e-12;
end
if isempty(RB_FAIL)
    RB_FAIL = 0;
end
if isempty(RB_TOTAL)
    RB_TOTAL = 0;
end

RB_TOTAL = RB_TOTAL + 1;

if (ischar(expected) && strcmp(expected, 'error'))
    threw = false;
    try
        fhandle(varargin{:});
    catch
        threw = true;
    end
    pass = threw;
    report(testname, pass, 'expected an error');
    return
end

try
    if (ischar(expected) && strcmp(expected, 'noerror'))
        fhandle(varargin{:});
        pass = true;
        report(testname, pass, 'no error');
        return
    end
    out = fhandle(varargin{:});
    [pass, msg] = compare(out, expected, RB_RTOL, RB_ATOL);
    report(testname, pass, msg);
catch err
    pass = false;
    report(testname, pass, sprintf('threw: %s', err.message));
end

function [ok, why] = compare(actual, expected, rtol, atol)
ok = false;
why = '';
if isa(expected, 'function_handle')
    try
        ok = expected(actual);
        if ~ok
            why = 'predicate returned false';
        end
    catch err
        why = sprintf('predicate threw: %s', err.message);
    end
    return
end
if isnumeric(expected) || islogical(expected)
    if ~isequal(size(actual), size(expected))
        why = sprintf('size mismatch: got [%s], expected [%s]', ...
                      sprintf('%d ', size(actual)), sprintf('%d ', size(expected)));
        return
    end
    a = double(actual(:));
    e = double(expected(:));
    if any(isnan(a) ~= isnan(e))
        why = 'NaN pattern mismatch';
        return
    end
    a(isnan(a)) = 0;
    e(isnan(e)) = 0;
    diff = abs(a - e);
    tol  = atol + rtol .* max(abs(a), abs(e));
    bad  = find(diff > tol, 1);
    if isempty(bad)
        ok = true;
    else
        why = sprintf('numerical mismatch at idx %d: |%g-%g|=%g > %g', ...
                      bad, a(bad), e(bad), diff(bad), tol(max(numel(tol), 1)));
    end
    return
end
if ischar(expected)
    ok = ischar(actual) && strcmp(actual, expected);
    if ~ok
        why = 'string mismatch';
    end
    return
end
if iscell(expected)
    if ~iscell(actual) || ~isequal(size(actual), size(expected))
        why = 'cell shape mismatch';
        return
    end
    for k = 1:numel(expected)
        [ok2, why2] = compare(actual{k}, expected{k}, rtol, atol);
        if ~ok2
            why = sprintf('cell{%d}: %s', k, why2);
            return
        end
    end
    ok = true;
    return
end
if isstruct(expected)
    if ~isstruct(actual)
        why = 'expected struct';
        return
    end
    fn = fieldnames(expected);
    for k = 1:numel(fn)
        if ~isfield(actual, fn{k})
            why = sprintf('missing field %s', fn{k});
            return
        end
        [ok2, why2] = compare(actual.(fn{k}), expected.(fn{k}), rtol, atol);
        if ~ok2
            why = sprintf('field %s: %s', fn{k}, why2);
            return
        end
    end
    ok = true;
    return
end
why = sprintf('unsupported expected type %s', class(expected));

function report(testname, pass, msg)
global RB_FAIL
if pass
    fprintf(1, 'Testing %s: ok\n', testname);
else
    RB_FAIL = RB_FAIL + 1;
    fprintf(2, 'Testing %s: FAIL  (%s)\n', testname, msg);
end
