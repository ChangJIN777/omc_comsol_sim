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
                % An unrecognised name is an ERROR, never a silent fall-through
                % to auto-detection: a typo would otherwise quietly write to a
                % different store than the one asked for, and the mistake would
                % only surface as "my resume found nothing".
                backend = char(backend);
                known = {'dbtoolbox', 'python', 'cli', 'jsonl'};
                if ~any(strcmp(backend, known))
                    error('OptimStore:unknownBackend', ...
                        ['Unknown backend "%s". Valid values are: %s.\n', ...
                         'Available on this machine: %s.'], ...
                        backend, strjoin(known, ', '), ...
                        strjoin(OptimStore.availableBackends(), ', '));
                end
                obj.backend = backend;
            end
            if strcmp(obj.backend, 'cli')
                obj.cliExe = OptimStore.findSqliteExe();
            end

            if obj.useJsonl()
                % No schema to create -- the JSONL file IS the store. Prove it
                % is writable now rather than at the first expensive solve.
                fid = fopen(obj.jsonPath, 'a');
                if fid < 0
                    error('OptimStore:jsonlUnwritable', ...
                        ['No SQLite backend is available and the JSONL ', ...
                         'fallback cannot be written either:\n  %s\n', ...
                         'Check the directory exists and is writable.'], ...
                        obj.jsonPath);
                end
                fclose(fid);
            else
                % A locked or networked filesystem can refuse SQLite outright.
                % Degrade to JSONL rather than abandoning the study.
                try
                    obj.initSchema();
                catch ME
                    warning('OptimStore:schemaFailed', ...
                        ['SQLite backend "%s" could not open %s (%s). ', ...
                         'Falling back to the JSONL store.'], ...
                        obj.backend, obj.dbPath, ME.message);
                    obj.backend = 'jsonl';
                end
            end
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
            if obj.useJsonl()
                obj.appendJsonl(struct('x_type', 'run', 'run_id', runId, ...
                    'created_ts', nowTs, 'updated_ts', nowTs, ...
                    'script', script, 'optimizer', optimizer, ...
                    'status', 'running', 'config_json', cfgJson));
                return;
            end
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
            if obj.useJsonl()
                % Append-only: the latest record for a run_id wins on read.
                obj.appendJsonl(struct('x_type', 'run', 'run_id', runId, ...
                    'updated_ts', OptimStore.nowEpoch(), 'status', status));
                return;
            end
            obj.exec('UPDATE runs SET status=?, updated_ts=? WHERE run_id=?', ...
                {status, OptimStore.nowEpoch(), runId});
        end

        function logEval(obj, runId, evalIdx, rec)
        %LOGEVAL  Persist one objective evaluation.
        %   rec fields (all optional except x): x = [dAR maxdef], oblong, wM,
        %   gOM, gMB, gPE, LSiV, Qmech, Qopt, lambdaNm, fitness, fitnessRaw,
        %   J, status, evLoc.  `fitness` is the gated score the Nelder-Mead
        %   ladder uses; `fitnessRaw` is the same figure of merit with the
        %   Q_mech gate DISABLED, which is what a GP surrogate needs -- a
        %   gated 0 tells it nothing about where the gate boundary lies.
        %   evLoc.  gMB/gPE are the moving-boundary and photoelastic parts
        %   of gOM; they are recorded at the same mode pair as gOM, so the
        %   three always correspond to one another.
        %   The JSONL mirror is written FIRST, so a SQLite failure cannot lose
        %   a solve that has already been paid for.
            x        = rec.x;
            paramKey = OptimStore.paramKeyOf(x);

            rec.run_id    = runId;
            rec.eval_idx  = evalIdx;
            rec.ts        = OptimStore.nowEpoch();
            rec.param_key = paramKey;   % so the JSONL reader can index on it
            rec.x_type    = 'eval';
            obj.appendJsonl(rec);

            if obj.useJsonl()
                return;     % the append above IS the write
            end

            g = @(f) OptimStore.getf(rec, f);
            try
                obj.exec(['INSERT OR REPLACE INTO evals (run_id, eval_idx, ts, ', ...
                          'param_key, dar, maxdef, oblong, wm_hz, gom_hz, ', ...
                          'gmb_hz, gpe_hz, ', ...
                          'lsiv_hz, q_mech, q_opt, lambda_nm, fitness, ', ...
                          'fitness_raw, ', ...
                          'objective_j, status, ev_loc) ', ...
                          'VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'], ...
                    {runId, evalIdx, rec.ts, paramKey, x(1), x(2), g('oblong'), ...
                     g('wM'), g('gOM'), g('gMB'), g('gPE'), ...
                     g('LSiV'), g('Qmech'), g('Qopt'), ...
                     g('lambdaNm'), g('fitness'), g('fitnessRaw'), ...
                     g('J'), g('status'), g('evLoc')});
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
            if obj.useJsonl()
                hit = obj.jsonlLookup(runId, x);
                return;
            end
            rows = obj.query(['SELECT eval_idx, dar, maxdef, oblong, fitness, ', ...
                              'objective_j, status, fitness_raw, q_mech ', ...
                              'FROM evals ', ...
                              'WHERE run_id=? AND param_key=? LIMIT 1'], ...
                {runId, OptimStore.paramKeyOf(x)});
            if isempty(rows)
                hit = [];
                return;
            end
            % fitnessRaw and Qmech are returned so that an optimizer with a
            % different scoring rule than the one that wrote the row can
            % still re-derive its own objective from the cached physics.
            hit = struct('eval_idx',   OptimStore.num(rows{1,1}), ...
                         'x',          [OptimStore.num(rows{1,2}), ...
                                        OptimStore.num(rows{1,3})], ...
                         'oblong',     OptimStore.num(rows{1,4}), ...
                         'fitness',    OptimStore.num(rows{1,5}), ...
                         'J',          OptimStore.num(rows{1,6}), ...
                         'status',     OptimStore.str(rows{1,7}), ...
                         'fitnessRaw', OptimStore.num(rows{1,8}), ...
                         'Qmech',      OptimStore.num(rows{1,9}));
        end

        function T = loadEvals(obj, runId)
        %LOADEVALS  All evaluations for a run, oldest first, as a table.
            if obj.useJsonl()
                T = obj.jsonlLoadEvals(runId);
                return;
            end
            rows = obj.query(['SELECT eval_idx, ts, dar, maxdef, oblong, ', ...
                              'wm_hz, gom_hz, gmb_hz, gpe_hz, ', ...
                              'lsiv_hz, q_mech, q_opt, ', ...
                              'lambda_nm, fitness, fitness_raw, ', ...
                              'objective_j, status ', ...
                              'FROM evals WHERE run_id=? ORDER BY eval_idx'], ...
                {runId});
            vars = {'eval','ts','dAR','maxdef','oblong','wM','gOM', ...
                    'gMB','gPE','LSiV', ...
                    'Qmech','Qopt','lambdaNm','fitness','fitnessRaw', ...
                    'J','status'};
            T = OptimStore.rowsToTable(rows, vars, numel(vars) - 1);
        end

        function [x, fitness, evalIdx] = best(obj, runId)
        %BEST  Highest-fitness evaluation of a run (x = [] if none succeeded).
            if obj.useJsonl()
                x = [];  fitness = NaN;  evalIdx = NaN;
                T = obj.jsonlLoadEvals(runId);
                if height(T) == 0
                    return;
                end
                ok = isfinite(T.fitness) & T.fitness > 0;
                if ~any(ok)
                    return;
                end
                idx = find(ok);
                [fitness, k] = max(T.fitness(ok));
                k = idx(k);
                x       = [T.dAR(k), T.maxdef(k)];
                evalIdx = T.eval(k);
                return;
            end
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
            if obj.useJsonl()
                n = height(obj.jsonlLoadEvals(runId));
                return;
            end
            rows = obj.query('SELECT COUNT(*) FROM evals WHERE run_id=?', {runId});
            n = OptimStore.num(rows{1,1});
        end

        function T = listRuns(obj)
        %LISTRUNS  Every run in the database, most recent activity first.
            if obj.useJsonl()
                T = obj.jsonlListRuns();
                return;
            end
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
                      'wm_hz REAL, gom_hz REAL, gmb_hz REAL, gpe_hz REAL, ', ...
                      'lsiv_hz REAL, q_mech REAL, ', ...
                      'q_opt REAL, lambda_nm REAL, fitness REAL, ', ...
                      'fitness_raw REAL, ', ...
                      'objective_j REAL, status TEXT, ev_loc TEXT, ', ...
                      'PRIMARY KEY (run_id, eval_idx))'], {});
            obj.exec(['CREATE INDEX IF NOT EXISTS idx_evals_lookup ', ...
                      'ON evals(run_id, param_key)'], {});
            obj.migrateSchema();
        end

        function migrateSchema(obj)
        %MIGRATESCHEMA  Add columns a database from an older version lacks.
        %   SQLite has no 'ADD COLUMN IF NOT EXISTS', and PRAGMA table_info
        %   is not uniformly readable across all three backends, so the
        %   ALTER is simply attempted and a duplicate-column error treated
        %   as success.  ALTER TABLE ADD COLUMN rewrites no rows, so this is
        %   cheap and non-destructive: pre-existing evaluations keep their
        %   data and get NULL for the new columns.
            newCols = {'gmb_hz REAL', 'gpe_hz REAL', 'fitness_raw REAL'};
            for k = 1:numel(newCols)
                try
                    obj.exec(['ALTER TABLE evals ADD COLUMN ', ...
                        newCols{k}], {});
                catch ME
                    if ~contains(lower(ME.message), 'duplicate column')
                        warning('OptimStore:migrate', ...
                            'Could not add column "%s": %s', ...
                            newCols{k}, ME.message);
                    end
                end
            end
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

        function tf = useJsonl(obj)
        %USEJSONL  True when the dependency-free JSONL store is in use.
            tf = strcmp(obj.backend, 'jsonl');
        end

        function [evals, runs] = jsonlRead(obj)
        %JSONLREAD  Parse the append-only mirror into eval and run records.
        %   Append-only means a key can appear more than once; the LAST record
        %   wins, which reproduces SQLite's INSERT OR REPLACE / UPDATE. Run
        %   records are merged field-by-field so a bare finishRun update does
        %   not erase the script/optimizer written by startRun.
            evals = {};
            runs  = struct();
            if exist(obj.jsonPath, 'file') ~= 2
                return;
            end

            txt = fileread(obj.jsonPath);
            lines = regexp(txt, '\r?\n', 'split');

            evalKeys = {};
            for i = 1:numel(lines)
                ln = strtrim(lines{i});
                if isempty(ln)
                    continue;
                end
                try
                    r = jsondecode(ln);
                catch
                    % A torn final line (killed mid-write) must not sink the
                    % whole history -- skip it and keep the rest.
                    continue;
                end
                if ~isstruct(r)
                    continue;
                end

                % Records written before x_type existed are all evaluations.
                kind = 'eval';
                if isfield(r, 'x_type')
                    kind = OptimStore.str(r.x_type);
                end

                switch kind
                    case 'run'
                        if ~isfield(r, 'run_id'); continue; end
                        f = matlab.lang.makeValidName(OptimStore.str(r.run_id));
                        if isfield(runs, f)
                            runs.(f) = OptimStore.mergeStruct(runs.(f), r);
                        else
                            runs.(f) = r;
                        end
                    otherwise
                        if ~isfield(r, 'run_id') || ~isfield(r, 'eval_idx')
                            continue;
                        end
                        key = sprintf('%s|%d', OptimStore.str(r.run_id), ...
                            round(OptimStore.num(r.eval_idx)));
                        prev = find(strcmp(evalKeys, key), 1);
                        if isempty(prev)
                            evalKeys{end+1} = key;   %#ok<AGROW>
                            evals{end+1}    = r;     %#ok<AGROW>
                        else
                            evals{prev} = r;         % last write wins
                        end
                end
            end
        end

        function hit = jsonlLookup(obj, runId, x)
        %JSONLLOOKUP  lookup() against the JSONL store.
            hit = [];
            key = OptimStore.paramKeyOf(x);
            evals = obj.jsonlRead();
            for i = 1:numel(evals)
                r = evals{i};
                if ~strcmp(OptimStore.str(r.run_id), runId)
                    continue;
                end
                % Prefer the stored key; recompute it for legacy rows.
                if isfield(r, 'param_key')
                    rk = OptimStore.str(r.param_key);
                elseif isfield(r, 'x')
                    rk = OptimStore.paramKeyOf(r.x);
                else
                    continue;
                end
                if ~strcmp(rk, key)
                    continue;
                end
                g = @(f) OptimStore.getf(r, f);
                xv = [NaN NaN];
                if isfield(r, 'x') && numel(r.x) >= 2
                    xv = [OptimStore.num(r.x(1)), OptimStore.num(r.x(2))];
                end
                hit = struct('eval_idx',   OptimStore.num(r.eval_idx), ...
                             'x',          xv, ...
                             'oblong',     OptimStore.num(g('oblong')), ...
                             'fitness',    OptimStore.num(g('fitness')), ...
                             'J',          OptimStore.num(g('J')), ...
                             'status',     OptimStore.str(g('status')), ...
                             'fitnessRaw', OptimStore.num(g('fitnessRaw')), ...
                             'Qmech',      OptimStore.num(g('Qmech')));
                return;
            end
        end

        function T = jsonlLoadEvals(obj, runId)
        %JSONLLOADEVALS  loadEvals() against the JSONL store.
            vars = {'eval','ts','dAR','maxdef','oblong','wM','gOM', ...
                    'gMB','gPE','LSiV', ...
                    'Qmech','Qopt','lambdaNm','fitness','fitnessRaw', ...
                    'J','status'};
            evals = obj.jsonlRead();
            keep  = false(1, numel(evals));
            for i = 1:numel(evals)
                keep(i) = strcmp(OptimStore.str(evals{i}.run_id), runId);
            end
            evals = evals(keep);
            if isempty(evals)
                T = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
                return;
            end

            n = numel(evals);
            C = cell(n, numel(vars));
            for i = 1:n
                r = evals{i};
                g = @(f) OptimStore.num(OptimStore.getf(r, f));
                xv = [NaN NaN];
                if isfield(r, 'x') && numel(r.x) >= 2
                    xv = [OptimStore.num(r.x(1)), OptimStore.num(r.x(2))];
                end
                C(i, :) = {OptimStore.num(r.eval_idx), g('ts'), ...
                           xv(1), xv(2), g('oblong'), g('wM'), g('gOM'), ...
                           g('gMB'), g('gPE'), g('LSiV'), ...
                           g('Qmech'), g('Qopt'), g('lambdaNm'), ...
                           g('fitness'), g('fitnessRaw'), g('J'), ...
                           OptimStore.str(OptimStore.getf(r, 'status'))};
            end
            T = cell2table(C, 'VariableNames', vars);
            for k = 1:(numel(vars) - 1)
                if iscell(T.(vars{k}))
                    T.(vars{k}) = cell2mat(T.(vars{k}));
                end
            end
            T = sortrows(T, 'eval');
        end

        function T = jsonlListRuns(obj)
        %JSONLLISTRUNS  listRuns() against the JSONL store.
            vars = {'run_id','script','optimizer','status','created', ...
                    'updated','nEvals'};
            [evals, runs] = obj.jsonlRead();
            f = fieldnames(runs);
            if isempty(f)
                T = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
                return;
            end
            C = cell(numel(f), numel(vars));
            for i = 1:numel(f)
                r  = runs.(f{i});
                id = OptimStore.str(r.run_id);
                nE = 0;
                for j = 1:numel(evals)
                    nE = nE + strcmp(OptimStore.str(evals{j}.run_id), id);
                end
                C(i, :) = {id, ...
                           OptimStore.str(OptimStore.getf(r, 'script')), ...
                           OptimStore.str(OptimStore.getf(r, 'optimizer')), ...
                           OptimStore.str(OptimStore.getf(r, 'status')), ...
                           OptimStore.num(OptimStore.getf(r, 'created_ts')), ...
                           OptimStore.num(OptimStore.getf(r, 'updated_ts')), ...
                           nE};
            end
            T = cell2table(C, 'VariableNames', vars);
            for k = 5:7
                if iscell(T.(vars{k}))
                    T.(vars{k}) = cell2mat(T.(vars{k}));
                end
            end
            T = sortrows(T, 'updated', 'descend');
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
            % Make NULLs visible: printed as '' they are indistinguishable
            % from a trailing separator that got trimmed away, which silently
            % shortens the row. See parseCliRows.
            fprintf(fid, '.nullvalue %s\n', OptimStore.cliNullToken());
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
        %DETECTBACKEND  Pick the best available storage backend on this machine.
        %   Tiers in order of preference, ending in one that cannot fail.
        %   This NEVER errors: a store that refuses to open would abandon a
        %   study that may already be COMSOL-hours deep, which is the exact
        %   outcome the whole class exists to prevent.
            if exist('sqlite', 'file') == 2 || exist('sqlite', 'builtin') == 5
                backend = 'dbtoolbox';
                return;
            end
            if OptimStore.probePython()
                backend = 'python';
                return;
            end
            if ~isempty(OptimStore.findSqliteExe())
                backend = 'cli';
                return;
            end
            % Last resort: the append-only JSONL mirror, read back in memory.
            % Base MATLAB only, no toolbox, no interpreter, no executable --
            % so it works on any machine the optimizer itself runs on.
            backend = 'jsonl';
            warning('OptimStore:jsonlFallback', ...
                ['No SQLite backend found on this machine, falling back to ', ...
                 'the dependency-free JSONL store.\n', ...
                 '  Database Toolbox sqlite() : %d\n', ...
                 '  MATLAB Python interface   : %s\n', ...
                 '  sqlite3 executable        : %s\n', ...
                 'Runs, resume and seeding all work; only ad-hoc SQL ', ...
                 'querying of the results is unavailable. To get SQLite ', ...
                 'back, configure pyenv or put sqlite3 on the PATH.'], ...
                exist('sqlite', 'file') == 2, ...
                OptimStore.pyStatusText(), 'not found');
        end

        function s = pyStatusText()
        %PYSTATUSTEXT  One-line description of the MATLAB Python interface.
        %   Warnings are silenced: querying pyenv on a machine with a broken
        %   or terminated interpreter can emit one, and this function exists
        %   only to DESCRIBE the situation, never to add noise to it.
            ws = warning('off', 'all');
            restoreW = onCleanup(@() warning(ws));
            try
                pe = pyenv;
                if strcmp(char(pe.Version), '')
                    s = 'not configured';
                else
                    s = sprintf('%s (%s)', char(pe.Version), char(pe.Status));
                end
            catch
                s = 'unavailable';
            end
        end

        function exe = findSqliteExe()
        %FINDSQLITEEXE  Locate a sqlite3 executable; '' if there is none.
        %   Searches the PATH first, then the usual install locations on each
        %   platform. MATLAB's system() does not always inherit the same PATH
        %   as an interactive shell, so a sqlite3 you can run in a terminal is
        %   not necessarily one system() can find -- hence the explicit list.
            exe = '';
            [st, ~] = system('sqlite3 -version');
            if st == 0
                exe = 'sqlite3';
                return;
            end

            home = getenv('USERPROFILE');       % Windows
            if isempty(home)
                home = getenv('HOME');          % Linux / macOS
            end

            if ispc
                cands = { ...
                    fullfile(home, 'anaconda3', 'Library', 'bin', 'sqlite3.exe'), ...
                    fullfile(home, 'miniconda3', 'Library', 'bin', 'sqlite3.exe'), ...
                    fullfile(getenv('CONDA_PREFIX'), 'Library', 'bin', 'sqlite3.exe'), ...
                    fullfile(getenv('ProgramFiles'), 'sqlite3', 'sqlite3.exe'), ...
                    fullfile(getenv('SystemRoot'), 'System32', 'sqlite3.exe')};
            else
                cands = { ...
                    '/usr/bin/sqlite3', ...
                    '/usr/local/bin/sqlite3', ...
                    '/opt/homebrew/bin/sqlite3', ...
                    '/opt/local/bin/sqlite3', ...
                    fullfile(home, 'anaconda3', 'bin', 'sqlite3'), ...
                    fullfile(home, 'miniconda3', 'bin', 'sqlite3'), ...
                    fullfile(getenv('CONDA_PREFIX'), 'bin', 'sqlite3')};
            end

            for k = 1:numel(cands)
                c = cands{k};
                % A missing env var collapses a candidate to a relative stub;
                % skip those rather than probing the working directory.
                if isempty(c) || exist(c, 'file') ~= 2
                    continue;
                end
                exe = c;
                return;
            end
        end

        function key = paramKeyOf(x)
        %PARAMKEYOF  Deterministic lookup key for a design point.
        %   Rounded to 10 significant digits so a value that has round-tripped
        %   through SQLite's REAL storage still matches on the way back in.
            key = sprintf('%.10g|%.10g', x(1), x(2));
        end

        function backends = availableBackends()
        %AVAILABLEBACKENDS  Every storage backend usable on this machine.
        %   Always contains at least 'jsonl', which needs nothing but a
        %   writable directory. Use it to see what a given workstation has.
            backends = {};
            if exist('sqlite', 'file') == 2 || exist('sqlite', 'builtin') == 5
                backends{end+1} = 'dbtoolbox';
            end
            if OptimStore.probePython()
                backends{end+1} = 'python';
            end
            if ~isempty(OptimStore.findSqliteExe())
                backends{end+1} = 'cli';
            end
            backends{end+1} = 'jsonl';
        end

        function ok = probePython()
        %PROBEPYTHON  Can we actually DO SQLite through the Python interface?
        %
        %   Importing sqlite3 is not a sufficient test. MATLAB's Python bridge
        %   can import a module and then fail on the next call -- most often
        %   with an out-of-process interpreter that has died, which surfaces as
        %   "Unable to communicate with the Python interpreter". Detecting on
        %   the import alone therefore advertises a backend that breaks later,
        %   mid-study, after real solve time has been spent.
        %
        %   So probe with a complete round-trip against an in-memory database.
        %   Warnings are silenced as well as errors: a terminating Python
        %   process WARNS rather than throwing, so a bare try/catch lets the
        %   message through and it looks like a failure even when we recover.
        %   The result is cached for the MATLAB session. On a machine where the
        %   interpreter cannot even be launched -- the Microsoft Store build of
        %   Python is the common case, since the ACLs on
        %   %ProgramFiles%\WindowsApps deny process creation -- each attempt
        %   costs a failed process spawn and prints OS-level noise. Probing
        %   once per session keeps that to a single occurrence instead of one
        %   per store construction.
            persistent cached
            if ~isempty(cached)
                ok = cached;
                return;
            end

            ok = false;
            ws = warning('off', 'all');
            restoreW = onCleanup(@() warning(ws));
            try
                conn = py.sqlite3.connect(':memory:');
                conn.execute('CREATE TABLE probe (a INTEGER, b TEXT)');
                conn.execute('INSERT INTO probe VALUES (?,?)', ...
                    py.tuple({py.int(int64(1)), py.str('x')}));
                cur  = conn.execute('SELECT a, b FROM probe');
                rows = cell(cur.fetchall());
                conn.close();
                ok = isscalar(rows);
            catch
                % Any failure at all -- not configured, wrong version, dead
                % interpreter, blocked by policy -- means "not usable here".
                ok = false;
            end
            cached = ok;
        end

        function ok = selfTestAll()
        %SELFTESTALL  Run selfTest against EVERY backend this machine has.
        %   The one-call answer to "will this work here?". Reports what is
        %   available and exercises each, so a machine missing SQLite entirely
        %   still gets its JSONL fallback verified.
            backends = OptimStore.availableBackends();
            fprintf('Backends available here: %s\n', strjoin(backends, ', '));
            fprintf('  Database Toolbox sqlite() : %d\n', ...
                exist('sqlite', 'file') == 2 || exist('sqlite', 'builtin') == 5);
            fprintf('  MATLAB Python interface   : %s\n', OptimStore.pyStatusText());
            exe = OptimStore.findSqliteExe();
            if isempty(exe); exe = 'not found'; end
            fprintf('  sqlite3 executable        : %s\n\n', exe);

            ok     = true;
            failed = {};
            for i = 1:numel(backends)
                try
                    OptimStore.selfTest([], backends{i});
                catch ME
                    % One broken backend must not stop the others: the point
                    % of this call is to find out what DOES work here.
                    ok = false;
                    failed{end+1} = backends{i}; %#ok<AGROW>
                    fprintf('  backend "%s" FAILED: %s\n', ...
                        backends{i}, ME.message);
                end
            end

            if ok
                fprintf('\nALL BACKENDS PASSED\n');
            else
                fprintf(['\n%d backend(s) failed: %s\n', ...
                         'This is not fatal -- set OPT.storeBackend to one ', ...
                         'that passed. "jsonl" needs nothing but a writable ', ...
                         'directory and is always available.\n'], ...
                    numel(failed), strjoin(failed, ', '));
            end
        end

        function ok = selfTest(dbPath, backend)
        %SELFTEST  Exercise a backend end-to-end. No COMSOL needed.
        %   ok = OptimStore.selfTest                 auto-detected backend
        %   ok = OptimStore.selfTest(dbPath)         a database you name
        %   ok = OptimStore.selfTest([], backend)    force one backend
        %
        %   See also OptimStore.selfTestAll, which runs every backend the
        %   machine supports.
            if nargin < 2
                backend = '';
            end
            if nargin < 1 || isempty(dbPath)
                tag = char(datetime('now', 'Format', 'yyyyMMddHHmmssSSS'));
                if ~isempty(backend)
                    tag = [tag, '_', backend];
                end
                dbPath = fullfile(tempdir, ...
                    ['optimstore_selftest_', tag, '.sqlite3']);
            end
            ok = false;
            store = OptimStore(dbPath, backend);
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
            % Reopen with the SAME backend: auto-detection would pick a
            % different store than the one just written to, and "reopen lost
            % everything" would then be a test artifact, not a real failure.
            store2 = OptimStore(dbPath, store.backend);
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

        function a = mergeStruct(a, b)
        %MERGESTRUCT  Overlay b's fields onto a, so a partial append-only
        %   update (e.g. finishRun writing only status/updated_ts) does not
        %   drop fields an earlier record supplied.
            f = fieldnames(b);
            for i = 1:numel(f)
                a.(f{i}) = b.(f{i});
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
        %
        %   Two traps here, both of which silently TRUNCATE a row rather than
        %   erroring, so a SELECT whose trailing columns are NULL comes back
        %   short and the caller indexes past the end:
        %
        %   1. strtrim removes every character with code <= 32 -- which
        %      includes char(31), the column separator. Trimming a line that
        %      ends in separators (the signature of trailing NULL columns)
        %      therefore deletes them. Only \r and \n may be stripped here.
        %   2. sqlite3 prints NULL as an empty string by default, so those
        %      trailing columns are empty AND their separators are the only
        %      evidence they existed. runCli sets .nullvalue to a sentinel so
        %      a NULL is a visible token; it is mapped back to '' below.
            raw = regexprep(raw, '\r', '');
            raw = regexprep(raw, '\n+$', '');    % trailing newlines only
            if isempty(raw)
                rows = {};
                return;
            end
            lines = strsplit(raw, newline);
            rows  = {};
            r     = 0;
            for i = 1:numel(lines)
                if isempty(lines{i})
                    continue;
                end
                parts = strsplit(lines{i}, char(31), 'CollapseDelimiters', false);
                r = r + 1;
                for k = 1:numel(parts)
                    if strcmp(parts{k}, OptimStore.cliNullToken())
                        rows{r, k} = ''; %#ok<AGROW>
                    else
                        rows{r, k} = parts{k}; %#ok<AGROW>
                    end
                end
            end
        end

        function t = cliNullToken()
        %CLINULLTOKEN  Stand-in the sqlite3 CLI prints for a NULL value.
        %   Must be something no real stored value can be; the status strings
        %   and paths this class stores are all plain text.
            t = '<<NULL>>';
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
