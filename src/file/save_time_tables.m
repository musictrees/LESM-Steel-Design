function save_time_tables(model,path_name)
    %--------------------------------------------------------------------------
    % Dynamic time functions (load variation patterns over time)
    include_constants;
    ntfcn = 0;
    timeFcns = struct('fcn',[],'next',[],'prev',[]);
    for n = 1:model.nnp
        % Check if node has dynamic loads
        fptr = model.nodes(n).load.getFcn();
        while ~isempty(fptr)
            if fptr.type == TABLE                
                % Check if time fcn is first on local chained list of structs
                if isempty(timeFcns.fcn)
                    ntfcn = ntfcn + 1;
                    timeFcns.fcn = fptr;
                else
                    % Go through current list of fcns to check if this fcn has
                    % already been read
                    ptr = timeFcns;
                    fcnAlreadyRead = false;
                    while ~isempty(ptr)
                        if ptr.fcn == fptr
                            fcnAlreadyRead = true;
                            break
                        end
                        ptr = ptr.prev;
                    end
                    if ~fcnAlreadyRead
                        ntfcn = ntfcn + 1;
                        timeFcns = struct('fcn',fptr,'next',[],'prev',timeFcns);
                        timeFcns.prev.next = timeFcns;
                    end
                end
            end
            fptr = fptr.next;
        end
    end
 ptr = timeFcns;
 for i=1:ntfcn
     save_time_table(ptr.fcn,path_name);
     ptr = ptr.prev;
 end
end
function save_time_table(fcn,path_name)

if strcmp(fcn.src_file,'')
    t = clock;
    file_name = sprintf('%stime_table_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f.bin',...
                         path_name,t(1),t(2),t(3),t(4),t(5),t(6));
else
    [~,file_name,~] = fileparts(fcn.src_file);
    file_name = strcat(path_name,file_name,'.bin');
end
fcn.src_file = file_name;

fid = fopen(file_name,'w');
fwrite(fid,length(fcn.x),'double');
fwrite(fid,fcn.x,'double');
fwrite(fid,fcn.F,'double');
fclose(fid); 
end