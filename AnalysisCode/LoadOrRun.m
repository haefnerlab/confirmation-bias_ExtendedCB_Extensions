function varargout = LoadOrRun(func, args, save_file, varargin)
%LOADORRUN load precomputed results from a file, or compute and save
%them to the file.
%
% ... = LOADORRUN(func, {arg1, arg2, ..}, save_file, ..) looks for results
% in save_file (a path to a .mat file). If the file does not exist,
% computes func(arg1, arg2, ..) and saves those results to the file. Return
% values are identical to whatever func returns.
%
% ... = LOADORRUN(..., '-recompute') forces results to be recomputed.
%
% ... = LOADORRUN(..., '-verbose') prints useful information.
%
% S = LOADORRUN(..., '-struct', 'fieldxyz') (only valid if the func returns
% a single struct, S). Use this form when func adds a single field
% 'fieldxyz' to a struct. Only the absence of the field S.fieldxyz will
% trigger the evaluation of func.

%% Parse extra flags and remove them from varargin
flag_recompute = strcmpi('-recompute', varargin);
varargin = varargin(~flag_recompute);
recompute = any(flag_recompute);

flag_verbose = strcmpi('-verbose', varargin);
varargin = varargin(~flag_verbose);
verbose = any(flag_verbose);

flag_struct = strcmpi('-struct', varargin);
struct_mode = any(flag_struct);
if struct_mode
    % 'right shift' to get logical array that indexes the next argument
    flag_field = [false flag_struct(1:end-1)];
    fieldname = varargin{flag_field};
    % clear both '-struct' and the fieldname from varargin
    varargin = varargin(~(flag_struct | flag_field));
    
    if nargout > 1
        error('Only 1 return value is allowed with the ''-struct'' flag');
    end
end

if ~isempty(varargin), warning('LoadOrRun got unexpected args.'); end

results = cell(1, nargout);

%% Ensure that parent directory of save_file exists.
[save_dir, uid, ext] = fileparts(save_file);
% if the uid itself contains a '.' and does not end in '.mat', it will
% be split across 'uid' and 'ext'.
if ~strcmp(ext, '.mat'), uid = strcat(uid, ext); end

if ~isempty(save_dir) && ~exist(save_dir, 'dir')
    if verbose
        fprintf('Creating output directory: %s\n', save_dir);
    end
    mkdir(save_dir);
end

MAX_FILENAME_LENGTH = 255;
uid_short = maybe_hash(uid, MAX_FILENAME_LENGTH - 6);
save_file = fullfile(save_dir, [uid_short '.mat']);
id_save_file = fullfile(save_dir, [uid_short '.id.mat']);

%% Determine whether a call to func is needed.
do_compute = ~exist(save_file, 'file') || recompute;

if ~do_compute && struct_mode
    % If file exists but the output doesn't contain 'fieldxyz' (in 'struct
    % mode'), then we still need to set do_compute=true
    contents = load(save_file);
    if ~isfield(contents.varargout{1}, fieldname)
        do_compute = true;
        if verbose
            fprintf('Computing struct field ''%s''\n', fieldname);
        end
    end
end

% Check modification times.
if ~do_compute
    func_info = functions(func);
    file_info = dir(func_info.file);
    memo_info = dir(save_file);
    % If .m file changed more recently than the saved memo file, issue a
    % warning.
    if ~isempty(file_info) && file_info.datenum > memo_info.datenum
        warning(['Source file %s changed since the memo-file %s was last updated. ' ...
            'Delete the memo file if the output is affected!! FURTHER DELETE OTHER DEPENDENT ' ...
            'MEMO FILES AS THESE CANNOT ''SEE'' THE CHANGE IN THEIR DEPENDENCIES!!'], ...
            func_info.file, save_file);
    end
end

% Check for hash collision. Note that save_file might be large, so we
% separately save the full uid in the '.id.mat' file, which is very fast
% to load and verify.
if ~do_compute && exist(id_save_file, 'file')
    id_contents = load(id_save_file);
    if ~strcmp(id_contents.uid, uid)
        warning('Hash collision!! Original uids:\n\t%s\n\t%s', id_contents.uid, uid);
        do_compute = true;
    end
end

%% Call func or load precomputed results.
if do_compute
    % Call func(args) and capture as many return values as have been
    % requested by whoever called this function.
    if verbose
        fprintf('Calling %s with %d outputs...\t\n', func2str(func), nargout);
    end
    [results{:}] = func(args{:});
    
    % Sanity-check outputs when in struct mode
    if struct_mode
        if ~isstruct(results{1})
            error('func must return a struct when the ''-struct'' flag is used');
        elseif ~isfield(results{1}, fieldname)
            warning('In ''-struct'' mode, but func failed to add a field named ''%s''', fieldname);
        end
    end
    
    if verbose, fprintf('done. Saving.\n'); end
    % Save results to the file.
    save(save_file, 'results');
    save(id_save_file, 'uid', '-v7.3');
else
    if verbose, fprintf('Loading precomputed results...\t\n'); end
    contents = load(save_file);
    if verbose, fprintf('done.\n'); end
    results = contents.results;
end

varargout = results;
end

function [ uid ] = maybe_hash( string_uid, max_length )
%MAYBE_HASH hashes string_uid if its length is larger than max_length.
%max_length defaults to 250.

if nargin < 2, max_length = 250; end

if length(string_uid) > max_length
    % Hash the string.
    uid = sprintf('%X', string2hash(string_uid));
else
    % It's short enough; keep the string as-is.
    uid = string_uid;
end

end