function run(mode,tol)
% Constants
MODE_RUN  = uint8(1);
MODE_UPD  = uint8(2);
RESULT_OK = uint8(0);
FAIL_EXT  = uint8(1);
FAIL_OPEN = uint8(2);
FAIL_REF  = uint8(3);
FAIL_VER  = uint8(4);
FAIL_ANL  = uint8(5);
FAIL_COMP = uint8(6);

% Close running instances
if strcmp(class(findall(0)),'matlab.graphics.Graphics') %#ok<STISA>
    delete(findall(0));
end

% Print header
fprintf('==================================================================\n');
fprintf('              LESM - Linear Elements Structure Model              \n');
fprintf('                      Version 3.0 - May 2022                      \n');
fprintf('                         Regression Tests                         \n');
fprintf('==================================================================\n\n');

% Check input
if nargin < 2
    fprintf(2,'Not enough input arguments.\n');
    fprintf('\nExiting...\n');
    return;
end
if mode ~= MODE_RUN && mode ~= MODE_UPD
    fprintf(2,'Invalid mode.\n');
    fprintf('\nExiting...\n');
    return;
end

% Get input files
filter = {'*.lsm','LESM Model File (*.lsm)'};
title = 'LESM Testing - Input files';
[files,path] = uigetfile(filter,title,'MultiSelect','on');
if (isequal(files,0))
    fprintf(2,'No file selected.\n');
    fprintf('\nExiting...\n');
    return;
end
files = string(files);

% Send warning before updating reference results
if mode == MODE_UPD
    fprintf('This action will erase the reference results currently stored.\n');
    reply = input('Are you sure you want to continue (Y/N)?\n',"s");
    if ~strcmp(reply,'y')   &&...
       ~strcmp(reply,'Y')   &&...
       ~strcmp(reply,'yes') &&...
       ~strcmp(reply,'Yes') &&...
       ~strcmp(reply,'YES')
        fprintf('\nExiting...\n');
        return;
    else
        fprintf('\n');
    end
end

% Display selected files
n_files = length(files);
fprintf('Number of selected files: %d\n',n_files);
for i = 1:n_files
    fprintf('%d - %s\n',i,files(i));
end
fprintf('\n------------------------------------------------------------------\n\n');

% Vector of result flags
all_results = zeros(n_files,1);

% Run each file
tic
for i = 1:n_files
    % Clear variables for new analysis
    clearvars -except i files tol mode...
                      path n_files all_results...
                      MODE_RUN MODE_UPD...
                      RESULT_OK FAIL_EXT FAIL_OPEN FAIL_REF FAIL_VER FAIL_ANL FAIL_COMP
    include_constants;
    
    % Initialize result flag
    result = RESULT_OK;
    
    % Get current file name
    cur_file = files(i);
    
    % Check file extension
    if result == RESULT_OK
        [~,test_name,ext] = fileparts(cur_file);
        if ext ~= ".lsm"
            result = FAIL_EXT;
            all_results(i) = result;
        end
    end
    
    % Print test number and name
    fprintf('TEST NUMBER: %d\n',i);
    fprintf('TEST NAME:   %s\n',test_name);
    
    % Open test file
    if result == RESULT_OK
        fid_cur = fopen(cur_file,'rt');
        if fid_cur < 0
            result = FAIL_OPEN;
            all_results(i) = result;
        end
    end
    
    % Open reference results file
    if result == RESULT_OK
        ref_file = strcat(path,test_name,"_ref.txt");
        if mode == MODE_RUN
            fid_ref   = fopen(ref_file,'rt');
            file_info = dir(ref_file);
            if fid_ref < 0 || file_info.bytes == 0
                result = FAIL_REF;
                all_results(i) = result;
            end
        elseif mode == MODE_UPD
            fid_ref = fopen(ref_file,'wt');
            if fid_ref < 0
                result = FAIL_REF;
                all_results(i) = result;
            end
        end
    end
    
    % Check version
    if result == RESULT_OK
        model = Model();
        vs = readFile(fid_cur,model,0);
        if ~vs
            result = FAIL_VER;
            all_results(i) = result;
        end
    end
    
    % Perform analysis
    if result == RESULT_OK
        switch model.whichSolver
            case STATIC_LINEAR
                model.drv = Drv_LES(0,model);
            case DYNAMIC_NEWMARK_LINEAR
                model.drv = Drv_LED(DYNAMIC_NEWMARK_LINEAR,0,model);
            case DYNAMIC_MODALSUP_LINEAR
                model.drv = Drv_LED(DYNAMIC_MODALSUP_LINEAR,0,model);
        end
        if model.drv.process() ~= 1
            result = FAIL_ANL;
            all_results(i) = result;
        end
    end
    
    % Write reference results
    if result == RESULT_OK && mode == MODE_UPD
        if any(fopen('all'),fid_ref)
            format1 = sprintf('%%.%df ',tol);
            format2 = sprintf('%%.%df\n',tol);
            if model.whichSolver == STATIC_LINEAR
                if issparse(model.D)
                    cur_displ = full(model.D);
                else
                    cur_displ = model.D;
                end
                if issparse(model.F)
                    cur_force = full(model.F);
                else
                    cur_force = model.F;
                end
                fprintf(fid_ref,'%%NODAL_DISPLACEMENTS\n');
                fprintf(fid_ref,'%d\n',length(cur_displ));
                fprintf(fid_ref,format2,cur_displ);
                fprintf(fid_ref,'\n');
                fprintf(fid_ref,'%%NODAL_FORCES\n');
                fprintf(fid_ref,'%d\n',length(cur_force));
                fprintf(fid_ref,format2,cur_force);
            else
                cur_freq  = model.W;
                cur_modes = model.V;
                cur_displ = model.results.dynamicDispl(:,:,1);
                fprintf(fid_ref,'%%VIBRATION_FREQUENCIES\n');
                fprintf(fid_ref,'%d\n',length(cur_freq));
                fprintf(fid_ref,format2,cur_freq);
                fprintf(fid_ref,'\n');
                fprintf(fid_ref,'%%VIBRATION_MODES\n');
                fprintf(fid_ref,'%d\n',size(cur_modes,1));
                for j = 1:size(cur_modes,1)
                    fprintf(fid_ref,format1,cur_modes(j,:));
                    fprintf(fid_ref,'\n');
                end
                fprintf(fid_ref,'\n');
                fprintf(fid_ref,'%%NODAL_DISPLACEMENTS_TRANSIENT\n');
                fprintf(fid_ref,'%d\n',size(cur_displ,2));
                for j = 1:size(cur_displ,2)
                    fprintf(fid_ref,format1,cur_displ(:,j));
                    fprintf(fid_ref,'\n');
                end
            end
        else
            result = FAIL_REF;
            all_results(i) = result;
        end
    end
    
    % Read reference results
    ref_res_1 = [];
    ref_res_2 = [];
    ref_res_3 = [];
    if result == RESULT_OK && mode == MODE_RUN
        warning('off','MATLAB:deblank:NonStringInput');
        if any(fopen('all'),fid_ref)
            while ~feof(fid_ref)
                line = deblank(fgetl(fid_ref));
                switch line
                    case '%NODAL_DISPLACEMENTS'
                        if issparse(model.D)
                            ref_res_1 = zeros(size(full(model.D)));
                        else
                            ref_res_1 = zeros(size(model.D));
                        end
                        num_lines = fscanf(fid_ref,'%d',1);
                        for j = 1:num_lines
                            ref_res_1(j,1) = fscanf(fid_ref,'%f',1);
                        end
                    case '%NODAL_FORCES'
                        if issparse(model.F)
                            ref_res_2 = zeros(size(full(model.F)));
                        else
                            ref_res_2 = zeros(size(model.F));
                        end
                        num_lines = fscanf(fid_ref,'%d',1);
                        for j = 1:num_lines
                            ref_res_2(j,1) = fscanf(fid_ref,'%f',1);
                        end
                    case '%VIBRATION_FREQUENCIES'
                        cur_freq  = model.W;
                        ref_res_1 = zeros(size(cur_freq));
                        num_lines = fscanf(fid_ref,'%d',1);
                        for j = 1:num_lines
                            ref_res_1(j,1) = fscanf(fid_ref,'%f',1);
                        end
                    case '%VIBRATION_MODES'
                        cur_modes = model.V;
                        ref_res_2 = zeros(size(cur_modes));
                        num_lines = fscanf(fid_ref,'%d',1);
                        for j = 1:num_lines
                            ref_res_2(j,:) = fscanf(fid_ref,'%f',size(cur_modes,2));
                        end
                    case '%NODAL_DISPLACEMENTS_TRANSIENT'
                        cur_displ = model.results.dynamicDispl(:,:,1);
                        ref_res_3 = zeros(size(cur_displ));
                        num_lines = fscanf(fid_ref,'%d',1);
                        for j = 1:num_lines
                            ref_res_3(:,j) = fscanf(fid_ref,'%f',size(cur_displ,1));
                        end
                end
            end
        else
            result = FAIL_REF;
            all_results(i) = result;
        end
    end
    
    % Compare results
    if result == RESULT_OK && mode == MODE_RUN
        if model.whichSolver == STATIC_LINEAR
            cur_displ = model.D;
            cur_force = model.F;
            if ~isequal(size(cur_displ),size(ref_res_1)) ||...
               ~all(abs(cur_displ-ref_res_1) < tol,'all')
                result = FAIL_COMP;
                all_results(i) = result;
            end
            if ~isequal(size(cur_force),size(ref_res_2)) ||...
               ~all(abs(cur_force-ref_res_2) < tol,'all')
                result = FAIL_COMP;
                all_results(i) = result;
            end
        else
            cur_freq  = model.W;
            cur_modes = model.V;
            cur_displ = model.results.dynamicDispl(:,:,1);
            if ~isequal(size(cur_freq),size(ref_res_1))  ||...
               ~all(abs(cur_freq-ref_res_1)  < tol,'all')
                result = FAIL_COMP;
                all_results(i) = result;
            end
            if ~isequal(size(cur_modes),size(ref_res_2)) ||...
               ~all(abs(cur_modes-ref_res_2) < tol,'all')
                result = FAIL_COMP;
                all_results(i) = result;
            end
            if ~isequal(size(cur_displ),size(ref_res_3)) ||...
               ~all(abs(cur_displ-ref_res_3) < tol,'all')
                result = FAIL_COMP;
                all_results(i) = result;
            end
        end
    end
    
    % Close files
    fclose('all');
    
    % Print test result
    if result == RESULT_OK
        fprintf('RESULT: OK\n\n');
    elseif result == FAIL_EXT
        fprintf(2,'RESULT: Invalid file extension\n\n');
    elseif result == FAIL_OPEN
        fprintf(2,'RESULT: Error opening test file\n\n');
    elseif result == FAIL_REF
        fprintf(2,'RESULT: Error opening reference results file\n\n');
    elseif result == FAIL_VER
        fprintf(2,'RESULT: Incompatible file version\n\n');
    elseif result == FAIL_ANL
        fprintf(2,'RESULT: Analysis failed\n\n');
    elseif result == FAIL_COMP
        fprintf(2,'RESULT: Diferent results from reference\n\n');
    end
end

% Final message
fprintf('\nFinished!\n');
fprintf('Total time: %.3f s\n',toc);

if all(all_results == RESULT_OK)
    if mode == MODE_RUN
        fprintf('\nAll results are in agreement with references!\n');
    elseif mode == MODE_UPD
        fprintf('\nAll reference results were updated!\n');
    end
elseif mode == MODE_RUN
    fprintf('\nSome tests failed!\n');
    if any(all_results == FAIL_EXT)  ||...
       any(all_results == FAIL_OPEN) ||...
       any(all_results == FAIL_REF)  ||...
       any(all_results == FAIL_VER)
        fprintf('\nTests that could not be started or compared with reference:\n');
        fprintf('%d ',find(all_results==FAIL_EXT));
        fprintf('%d ',find(all_results==FAIL_OPEN));
        fprintf('%d ',find(all_results==FAIL_REF));
        fprintf('%d ',find(all_results==FAIL_VER));
        fprintf('\n');
    end
    if any(all_results == FAIL_ANL)
        fprintf('\nTests with failed analysis:\n');
        fprintf('%d ',find(all_results==FAIL_ANL));
        fprintf('\n');
    end
    if any(all_results == FAIL_COMP)
        fprintf('\nTests not matching reference results:\n');
        fprintf('%d ',find(all_results==FAIL_COMP));
        fprintf('\n');
    end
elseif mode == MODE_UPD
    fprintf('\nSome tests failed to update reference results:\n');
    if any(all_results == FAIL_EXT)  ||...
       any(all_results == FAIL_OPEN) ||...
       any(all_results == FAIL_REF)  ||...
       any(all_results == FAIL_VER)  ||...
       any(all_results == FAIL_ANL)
        fprintf('%d ',find(all_results==FAIL_EXT));
        fprintf('%d ',find(all_results==FAIL_OPEN));
        fprintf('%d ',find(all_results==FAIL_REF));
        fprintf('%d ',find(all_results==FAIL_VER));
        fprintf('%d ',find(all_results==FAIL_ANL));
        fprintf('\n');
    end
end
end