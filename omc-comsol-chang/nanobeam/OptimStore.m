classdef OptimStore < handle
%OPTIMSTORE  SQLite-backed persistence for stop/resume of long optimization runs.
%
%   Every objective evaluation in a COMSOL optimization costs minutes, so a run
%   that is interrupted -- a LiveLink crash, a reboot, an interactive Ctrl-C --
%   must not lose the solves it already paid for.  OptimStore writes one row per
%   evaluation to a SQLite database as it happens, so a later run can look up
%   any point that was already evaluated and skip the FEM call entirely.
%
%   The design mirrors python-scripts/src/database.py in this repo: one table
%   written immediately, plus an append-only JSONL mirror so a filesystem that
%   cannot do SQLite locking still never loses a result.
%
%   USAGE
%     store = OptimStore('./test/runs.sqlite3');
%     store.startRun('myrun_01', 'optimize_oblong_maxdef', 'neldermead', OPT);
%     store.logEval('myrun_01', 7, rec);              % rec: see logEval
%     hit = store.lookup('myrun_01', [1.2541 0.22]);  % [] if not evaluated
%     T   = store.loadEvals('myrun_01');              % everything so far
%     store.finishRun('myrun_01', 'complete');
%
%   BACKENDS (auto-detected, in this order)
%     'dbtoolbox' - MATLAB's built-in sqlite() object   (Database Toolbox)
%     'python'    - py.sqlite3 through MATLAB's Python interface
%     'cli'       - the sqlite3 executable driven via system()
%   Whichever is chosen, an append-only <db>.jsonl mirror is written too, so a
%   failed SQLite write never silently discards an expensive solve.
%
%   Verify the active backend end-to-end with:  OptimStore.selfTest

    properties (SetAccess = immutable)
        dbPath   char       % absolute path to the .sqlite3 file
        jsonPath char       % append-only mirror alongside it
        backend  char       % 'dbtoolbox' | 'python' | 'cli'
    end

    properties (Access = private)
        cliExe = 'sqlite3'  % resolved sqlite3 executable for the 'cli' backend
    end

    methods
        function obj = OptimStore(dbPath, backend)
        %OPTIMSTORE  Open (creating if needed) a run database.
        %   store = OptimStore(dbPath)           auto-detects the backend
        %   store = OptimStore(dbPath, backend)  forces one (mainly for tests)
            narginchk(1, 2);
            validateattributes(dbPath, {'char','string'}, {'scalartext'}, ...
                mfilename, 'dbPath');
            dbPath = char(dbPath);

            % Normalize separators the way RunNanobeamFEM does: a hard-coded
            % '\' is a literal character, not a separator, off Windows.
            dbPath = strrep(strrep(dbPath, '\', filesep), '/', filesep);
            dbDir  = fileparts(dbPath);
            if ~isempty(dbDir) && ~exist(dbDir, 'dir')
                mkdir(dbDir);
            end

            obj.dbPath   = dbPath;
            [p, n]       = fileparts(dbPath);
            obj.jsonPath = fullfile(p, [n, '.jsonl']);

            if nargin < 2 || isempty(backend)
                obj.backend = OptimStore.detectBackend();
            else
                obj.backend = char(backend);
            end
            if strcmp(obj.backend, 'cli')
                obj.cliExe = OptimStore.findSqliteExe();
            end

            obj.initSchema();
        end

        function startRun(obj, runId, script, optimizer, cfg)
        %STARTRUN  Register (or re-open) a run and snapshot its configuration.
        %   Re-running with an existing runId marks it 'running' again and
        %   leaves its evaluations in place -- that is what makes resume work.
            if nargin < 5
                cfg = struct();
            end
            cfgJson = '';
            try
                cfgJson = jsonencode(OptimStore.sanitizeForJson(cfg));
            catch ME
                warning('OptimStore:cfgEncode', ...
                    'Could not encode config to JSON: %s', ME.message);
            end
            nowTs = OptimStore.nowEpoch();
            obj.exec(['INSERT INTO runs (run_id, created_ts, updated_ts, ', ...
                      'script, optimizer, status, config_json) ', ...
                      'VALUES (?,?,?,?,?,?,?) ', ...
                      'ON CONFLICT(run_id) DO UPDATE SET ', ...
                      'updated_ts=excluded.updated_ts, status=''running'', ', ...
                      'optimizer=excluded.optimizer'], ...
                {runId, nowTs, nowTs, script, optimizer, 'running', cfgJson});
        end

        function finishRun(obj, runId, status)
        %FINISHRUN  Mark a run 'complete', 'interrupted', or 'failed'.
            obj.exec('UPDATE runs SET status=?, updated_ts=? WHERE run_id=?', ...
                {status, OptimStore.nowEpoch(), runId});
        end

        function logEval(obj, runId, evalIdx, rec)
        %LOGEVAL  Persist one objective evaluation.
        %   rec fields (all optional except x): x = [dAR maxdef], oblong, wM,
        %   gOM, LSiV, Qmech, Qopt, lambdaNm, fitness, J, status, evLoc.
        %   The JSONL mirror is written FIRST, so a SQLite failure cannot lose
        %   a solve that has already been paid for.
            x        = rec.x;
            paramKey = OptimStore.paramKeyOf(x);

            rec.run_id   = runId;
            rec.eval_idx = evalIdx;
            rec.ts       = OptimStore.nowEpoch();
            obj.appendJsonl(rec);

            g = @(f) OptimStore.getf(rec, f);
            try
                obj.exec(['INSERT OR REPLACE INTO evals (run_id, eval_idx, ts, ', ...
                          'param_key, dar, maxdef, oblong, wm_hz, gom_hz, ', ...
                          'lsiv_hz, q_mech, q_opt, lambda_nm, fitness, ', ...
                          'objective_j, status, ev_loc) ', ...
                          'VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'], ...
                    {runId, evalIdx, rec.ts, paramKey, x(1), x(2), g('oblong'), ...
                     g('wM'), g('gOM'), g('LSiV'), g('Qmech'), g('Qopt'), ...
                     g('lambdaNm'), g('fitness'), g('J'), g('status'), g('evLoc')});
                obj.exec('UPDATE runs SET updated_ts=? WHERE run_id=?', ...
                    {rec.ts, runId});
            catch ME
                % The JSONL mirror above already holds the row, so the solve is
                % not lost -- warn and continue rather than aborting the study.
                warning('OptimStore:writeFailed', ...
                    'SQLite write failed for eval %d (kept in %s): %s', ...
                    evalIdx, obj.jsonPath, ME.message);
            end
        end

        function hit = lookup(obj, runId, x)
        %LOOKUP  Return the stored evaluation at x, or [] if there is none.
        %   This is what makes resume cheap: a deterministic optimizer replays
        %   its path through already-evaluated points at zero FEM cost.
            rows = obj.query(['SELECT eval_idx, dar, maxdef, oblong, fitness, ', ...
                              'objective_j, status FROM evals ', ...
                              'WHERE run_id=? AND param_key=? LIMIT 1'], ...
                {runId, OptimStore.paramKeyOf(x)});
            if isempty(rows)
                hit = [];
                return;
            end
            hit = struct('eval_idx', OptimStore.num(rows{1,1}), ...
                         'x',        [OptimStore.num(rows{1,2}), ...
                                      OptimStore.num(rows{1,3})], ...
                         'oblong',   OptimStore.num(rows{1,4}), ...
                         'fitness',  OptimStore.num(rows{1,5}), ...
                         'J',        OptimStore.num(rows{1,6}), ...
                         'status',   OptimStore.str(rows{1,7}));
        end

        function T = loadEvals(obj, runId)
        %LOADEVALS  All evaluations for a run, oldest first, as a table.
            rows = obj.query(['SELECT eval_idx, ts, dar, maxdef, oblong, ', ...
                              'wm_hz, gom_hz, lsiv_hz, q_mech, q_opt, ', ...
                              'lambda_nm, fitness, objective_j, status ', ...
                              'FROM evals WHERE run_id=? ORDER BY eval_idx'], ...
                {runId});
            vars = {'eval','ts','dAR','maxdef','oblong','wM','gOM','LSiV', ...
                    'Qmech','Qopt','lambdaNm','fitness','J','status'};
            T = OptimStore.rowsToTable(rows, vars, numel(vars) - 1);
        end

        function [x, fitness, evalIdx] = best(obj, runId)
        %BEST  Highest-fitness evaluation of a run (x = [] if none succeeded).
            rows = obj.query(['SELECT dar, maxdef, fitness, eval_idx ', ...
                              'FROM evals WHERE run_id=? AND fitness IS NOT NULL ', ...
                              'AND fitness>0 ORDER BY fitness DESC LIMIT 1'], ...
                {runId});
            if isempty(rows)
                x = [];  fitness = NaN;  evalIdx = NaN;
                return;
            end
            x       = [OptimStore.num(rows{1,1}), OptimStore.num(rows{1,2})];
            fitness = OptimStore.num(rows{1,3});
            evalIdx = OptimStore.num(rows{1,4});
        end

        function n = evalCount(obj, runId)
        %EVALCOUNT  Number of evaluations already stored for a run.
            rows = obj.query('SELECT COUNT(*) FROM evals WHERE run_id=?', {runId});
            n = OptimStore.num(rows{1,1});
        end

        function T = listRuns(obj)
        %LISTRUNS  Every run in the database, most recent activity first.
            rows = obj.query(['SELECT run_id, script, optimizer, status, ', ...
                              'created_ts, updated_ts, ', ...
                              '(SELECT COUNT(*) FROM evals e WHERE ', ...
                              'e.run_id=r.run_id) FROM runs r ', ...
                              'ORDER BY updated_ts DESC'], {});
            vars = {'run_id','script','optimizer','status','created', ...
                    'updated','nEvals'};
            if isempty(rows)
                T = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
                return;
            end
            n = size(rows, 1);
            C = cell(n, numel(vars));
            for i = 1:n
                for k = 1:4
                    C{i, k} = OptimStore.str(rows{i, k});
                end
                for k = 5:7
                    C{i, k} = OptimStore.num(rows{i, k});
                end
            end
            T = cell2table(C, 'VariableNames', vars);
            for k = 5:7
                if iscell(T.(vars{k}))
                    T.(vars{k}) = cell2mat(T.(vars{k}));
                end
            end
        end
    end

    %% ------------------------------------------------------------------
    %  Backend plumbing
    %  ------------------------------------------------------------------
    methods (Access = private)
        function initSchema(obj)
            obj.exec(['CREATE TABLE IF NOT EXISTS runs (', ...
                      'run_id TEXT PRIMARY KEY, created_ts REAL, ', ...
                      'updated_ts REAL, script TEXT, optimizer TEXT, ', ...
                      'status TEXT, config_json TEXT)'], {});
            obj.exec(['CREATE TABLE IF NOT EXISTS evals (', ...
                      'run_id TEXT, eval_idx INTEGER, ts REAL, ', ...
                      'param_key TEXT, dar REAL, maxdef REAL, oblong REAL, ', ...
                      'wm_hz REAL, gom_hz REAL, lsiv_hz REAL, q_mech REAL, ', ...
                      'q_opt REAL, lambda_nm REAL, fitness REAL, ', ...
                      'objective_j REAL, status TEXT, ev_loc TEXT, ', ...
                      'PRIMARY KEY (run_id, eval_idx))'], {});
            obj.exec(['CREATE INDEX IF NOT EXISTS idx_evals_lookup ', ...
                      'ON evals(run_id, param_key)'], {});
        end

        function appendJsonl(obj, rec)
            fid = -1;
            try
                fid = fopen(obj.jsonPath, 'a');
                if fid < 0
                    return;
                end
                fprintf(fid, '%s\n', jsonencode(OptimStore.sanitizeForJson(rec)));
            catch
                % A mirror that cannot be written must never abort a study.
            end
            if fid >= 0
                fclose(fid);
            end
        end

        function exec(obj, sql, params)
        %EXEC  Run a statement that returns nothing.
            if nargin < 3
                params = {};
            end
            switch obj.backend
                case 'dbtoolbox'
                    conn = sqlite(obj.dbPath, 'connect');
                    cleanupObj = onCleanup(@() close(conn));
                    execute(conn, OptimStore.inlineParams(sql, params));
                case 'python'
                    conn = py.sqlite3.connect(obj.dbPath);
                    try
                        conn.execute(sql, OptimStore.toPyTuple(params));
                        conn.commit();
                    catch ME
                        OptimStore.tryClose(conn);
                        rethrow(ME);
                    end
                    OptimStore.tryClose(conn);
                case 'cli'
                    obj.runCli(OptimStore.inlineParams(sql, params), {});
                otherwise
                    error('OptimStore:backend', ...
                        'Unknown backend "%s"', obj.backend);
            end
        end

        function rows = query(obj, sql, params)
        %QUERY  Run a SELECT; returns an n-by-m cell array ({} if no rows).
            if nargin < 3
                params = {};
            end
            switch obj.backend
                case 'dbtoolbox'
                    conn = sqlite(obj.dbPath, 'connect');
                    cleanupObj = onCleanup(@() close(conn));
                    res = fetch(conn, OptimStore.inlineParams(sql, params));
                    if istable(res)
                        res = table2cell(res);
                    end
                    rows = res;
                case 'python'
                    conn = py.sqlite3.connect(obj.dbPath);
                    try
                        cur    = conn.execute(sql, OptimStore.toPyTuple(params));
                        pyRows = cell(cur.fetchall());
                    catch ME
                        OptimStore.tryClose(conn);
                        rethrow(ME);
                    end
                    OptimStore.tryClose(conn);
                    rows = {};
                    for i = 1:numel(pyRows)
                        r = cell(pyRows{i});
                        for k = 1:numel(r)
                            rows{i, k} = r{k}; %#ok<AGROW>
                        end
                    end
                case 'cli'
                    % A unit-separator (0x1F) column split is unambiguous even
                    % when a text field legitimately contains a comma or pipe.
                    raw = obj.runCli(OptimStore.inlineParams(sql, params), ...
                        {'-separator', char(31)});
                    rows = OptimStore.parseCliRows(raw);
                otherwise
                    error('OptimStore:backend', ...
                        'Unknown backend "%s"', obj.backend);
            end
        end

        function out = runCli(obj, sql, extraArgs)
            args = '';
            for i = 1:2:numel(extraArgs)
                args = [args, ' ', extraArgs{i}, ' "', extraArgs{i+1}, '"']; %#ok<AGROW>
            end
            % The statement goes in via stdin rather than on the command line,
            % so its quoting does not have to survive a second shell escaping.
            tmp = [tempname, '.sql'];
            fid = fopen(tmp, 'w');
            fprintf(fid, '%s;\n', sql);
            fclose(fid);
            cleanupObj = onCleanup(@() delete(tmp));
            cmd = sprintf('"%s"%s "%s" < "%s"', obj.cliExe, args, obj.dbPath, tmp);
            [st, out] = system(cmd);
            if st ~= 0
                error('OptimStore:cli', 'sqlite3 failed (%d): %s', st, strtrim(out));
            end
        end
    end

    %% ------------------------------------------------------------------
    %  Public statics
    %  ------------------------------------------------------------------
    methods (Static)
        function backend = detectBackend()
        %DETECTBACKEND  Pick the best available SQLite path on this machine.
            if exist('sqlite', 'file') == 2 || exist('sqlite', 'builtin') == 5
                backend = 'dbtoolbox';
                return;
            end
            try
                py.importlib.import_module('sqlite3');
                backend = 'python';
                return;
            catch
                % Python interface not configured -- fall through to the CLI.
            end
            if ~isempty(OptimStore.findSqliteExe())
                backend = 'cli';
                return;
            end
            error('OptimStore:noBackend', ...
                ['No SQLite backend available. Install Database Toolbox, ', ...
                 'configure MATLAB''s Python interface (see pyenv), or put ', ...
                 'sqlite3 on the system PATH.']);
        end

        function exe = findSqliteExe()
        %FINDSQLITEEXE  Locate a sqlite3 executable; '' if there is none.
            exe = '';
            [st, ~] = system('sqlite3 -version');
            if st == 0
                exe = 'sqlite3';
                return;
            end
            % Anaconda ships one, and is commonly present on these machines
            % without being on the PATH that MATLAB's system() inherits.
            cand = fullfile(getenv('USERPROFILE'), 'anaconda3', 'Library', ...
                'bin', 'sqlite3.exe');
            if exist(cand, 'file') == 2
                exe = cand;
            end
        end

        function key = paramKeyOf(x)
        %PARAMKEYOF  Deterministic lookup key for a design point.
        %   Rounded to 10 significant digits so a value that has round-tripped
        %   through SQLite's REAL storage still matches on the way back in.
            key = sprintf('%.10g|%.10g', x(1), x(2));
        end

        function ok = selfTest(dbPath)
        %SELFTEST  Exercise the active backend end-to-end. No COMSOL needed.
        %   ok = OptimStore.selfTest            uses a temporary database
        %   ok = OptimStore.selfTest(dbPath)    uses a database you name
            if nargin < 1
                dbPath = fullfile(tempdir, ['optimstore_selftest_', ...
                    datestr(now, 'yyyymmddHHMMSSFFF'), '.sqlite3']); %#ok<TNOW1,DATST>
            end
            ok = false;
            store = OptimStore(dbPath);
            fprintf('OptimStore.selfTest: backend = %s\n', store.backend);
            fprintf('  db  -> %s\n', store.dbPath);

            runId = 'selftest_run';
            store.startRun(runId, 'OptimStore.selfTest', 'neldermead', ...
                struct('note', 'synthetic', 'fitnessFcn', @(ds) ds));

            % Three synthetic evaluations: a good one, a gated one, a failure.
            store.logEval(runId, 1, struct('x', [1.1135, 0.22], ...
                'oblong', 1.96, 'Qmech', 2.5e7, 'Qopt', 5e5, ...
                'lambdaNm', 1551.2, 'fitness', 4.5e16, 'J', -16.653, ...
                'status', 'ok', 'evLoc', './eval_0001/'));
            store.logEval(runId, 2, struct('x', [1.2000, 0.23], ...
                'oblong', 1.80, 'Qmech', 3.1e6, 'Qopt', 4e5, ...
                'lambdaNm', 1560.0, 'fitness', 0, 'J', 1e3, ...
                'status', 'Qmech_gated', 'evLoc', './eval_0002/'));
            store.logEval(runId, 3, struct('x', [0.9000, 0.30], ...
                'fitness', NaN, 'J', 2e3, 'status', 'failed'));

            assert(store.evalCount(runId) == 3, 'evalCount wrong');

            hit = store.lookup(runId, [1.1135, 0.22]);
            assert(~isempty(hit), 'lookup missed a stored point');
            assert(abs(hit.J - (-16.653)) < 1e-9, 'lookup returned wrong J');
            assert(strcmp(hit.status, 'ok'), 'lookup returned wrong status');

            assert(isempty(store.lookup(runId, [0.5, 0.5])), ...
                'lookup invented a point that was never evaluated');

            [xb, fb] = store.best(runId);
            assert(abs(xb(1) - 1.1135) < 1e-9 && abs(fb - 4.5e16) < 1e6, ...
                'best() picked the wrong row');

            T = store.loadEvals(runId);
            assert(height(T) == 3, 'loadEvals row count wrong');
            assert(isnan(T.fitness(3)), 'NaN fitness did not round-trip as NULL');
            assert(strcmp(T.status{2}, 'Qmech_gated'), 'status did not round-trip');

            store.finishRun(runId, 'complete');
            R = store.listRuns();
            assert(any(strcmp(R.run_id, runId)), 'listRuns missed the run');
            assert(strcmp(R.status{strcmp(R.run_id, runId)}, 'complete'), ...
                'finishRun did not update status');

            % Reopening must see everything -- this is the actual resume path.
            store2 = OptimStore(dbPath);
            assert(store2.evalCount(runId) == 3, 'reopen lost evaluations');
            assert(~isempty(store2.lookup(runId, [1.2000, 0.23])), ...
                'reopen lost the lookup index');

            assert(exist(store.jsonPath, 'file') == 2, 'JSONL mirror missing');

            fprintf('  all assertions passed (3 evals, lookup, best, reopen)\n');
            ok = true;
        end
    end

    %% ------------------------------------------------------------------
    methods (Static, Access = private)
        function t = nowEpoch()
            t = posixtime(datetime('now', 'TimeZone', 'UTC'));
        end

        function tryClose(conn)
            try
                conn.close();
            catch
                % Already closed, or the interpreter is gone -- nothing to do.
            end
        end

        function v = getf(s, f)
            if isfield(s, f) && ~isempty(s.(f))
                v = s.(f);
            else
                v = NaN;
            end
        end

        function T = rowsToTable(rows, vars, nNumeric)
        %ROWSTOTABLE  Cell rows -> table, first nNumeric columns as doubles.
            if isempty(rows)
                T = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
                return;
            end
            n = size(rows, 1);
            C = cell(n, numel(vars));
            for i = 1:n
                for k = 1:nNumeric
                    C{i, k} = OptimStore.num(rows{i, k});
                end
                for k = (nNumeric + 1):numel(vars)
                    C{i, k} = OptimStore.str(rows{i, k});
                end
            end
            T = cell2table(C, 'VariableNames', vars);
            % cell2table already flattens a column whose entries are all
            % same-type scalars; only the ones it left as cells need help.
            for k = 1:nNumeric
                if iscell(T.(vars{k}))
                    T.(vars{k}) = cell2mat(T.(vars{k}));
                end
            end
        end

        function v = num(raw)
        %NUM  Coerce a backend scalar to double (NaN for SQL NULL / empty).
            if isa(raw, 'py.NoneType')
                v = NaN;
            elseif isnumeric(raw)
                if isempty(raw)
                    v = NaN;
                else
                    v = double(raw);
                end
            elseif isa(raw, 'py.int') || isa(raw, 'py.float')
                v = double(raw);
            else
                s = OptimStore.str(raw);
                if isempty(s)
                    v = NaN;
                else
                    v = str2double(s);
                end
            end
        end

        function s = str(raw)
        %STR  Coerce a backend scalar to char ('' for SQL NULL / empty).
            if isa(raw, 'py.NoneType')
                s = '';
            elseif ischar(raw)
                s = raw;
            elseif isstring(raw)
                s = char(raw);
            elseif isa(raw, 'py.str')
                s = char(raw);
            elseif isempty(raw)
                s = '';
            else
                s = char(string(raw));
            end
        end

        function t = toPyTuple(params)
        %TOPYTUPLE  MATLAB cell -> Python tuple, with NaN mapped to SQL NULL.
            c = cell(1, numel(params));
            for i = 1:numel(params)
                p = params{i};
                if isnumeric(p) && isscalar(p) && isnan(p)
                    c{i} = py.None;
                elseif isnumeric(p) && isscalar(p)
                    if p == fix(p) && abs(p) < 2^53
                        c{i} = py.int(int64(p));
                    else
                        c{i} = py.float(p);
                    end
                elseif ischar(p) || isstring(p)
                    c{i} = py.str(char(p));
                elseif isempty(p)
                    c{i} = py.None;
                else
                    c{i} = py.str(char(string(p)));
                end
            end
            t = py.tuple(c);
        end

        function sql = inlineParams(sql, params)
        %INLINEPARAMS  Substitute ? placeholders for the backends that have no
        %   native parameter binding (Database Toolbox execute/fetch, and the
        %   CLI).  Text is escaped by doubling single quotes; NaN becomes NULL.
            for i = 1:numel(params)
                p = params{i};
                if isnumeric(p) && isscalar(p) && isnan(p)
                    lit = 'NULL';
                elseif isnumeric(p) && isscalar(p)
                    lit = sprintf('%.17g', p);
                elseif isempty(p)
                    lit = 'NULL';
                else
                    lit = ['''', strrep(char(p), '''', ''''''), ''''];
                end
                k = strfind(sql, '?');
                if isempty(k)
                    error('OptimStore:paramCount', ...
                        'More parameters than ? placeholders in: %s', sql);
                end
                sql = [sql(1:k(1)-1), lit, sql(k(1)+1:end)];
            end
        end

        function rows = parseCliRows(raw)
        %PARSECLIROWS  Split sqlite3 CLI output into an n-by-m cell array.
            raw = regexprep(raw, '\r', '');
            if isempty(strtrim(raw))
                rows = {};
                return;
            end
            lines = strsplit(strtrim(raw), newline);
            rows  = {};
            r     = 0;
            for i = 1:numel(lines)
                if isempty(strtrim(lines{i}))
                    continue;
                end
                parts = strsplit(lines{i}, char(31), 'CollapseDelimiters', false);
                r = r + 1;
                for k = 1:numel(parts)
                    rows{r, k} = parts{k}; %#ok<AGROW>
                end
            end
        end

        function s = sanitizeForJson(s)
        %SANITIZEFORJSON  Replace fields jsonencode cannot represent (function
        %   handles, COMSOL Java objects) so a config snapshot never throws.
            if ~isstruct(s)
                return;
            end
            f = fieldnames(s);
            for i = 1:numel(f)
                v = s.(f{i});
                if isa(v, 'function_handle')
                    s.(f{i}) = func2str(v);
                elseif isjava(v)
                    s.(f{i}) = '<java object>';
                elseif isstruct(v)
                    s.(f{i}) = OptimStore.sanitizeForJson(v);
                end
            end
        end
    end
end
