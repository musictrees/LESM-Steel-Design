function [x,F] = read_time_table(file_name)
[~,~,ext] = fileparts(file_name);

if strcmp(ext,'.txt')
    
fid = fopen(file_name);
if fid < 0
    x = [];
    F = [];
    return
end

n = fscanf(fid,'%d',1);
x = zeros(n,1);
F = zeros(n,1);

cont = 1;

while ~feof(fid) && cont <= n
    v = fscanf(fid,'%f %f', [ 1 2]);
    x(cont) = v(1);
    F(cont) = v(2);
    cont = cont+1;
end

fclose(fid); 

elseif strcmp(ext,'.bin')
    fid = fopen(file_name,'r');
    
    if fid < 0
        x = [];
        F = [];
        return
    end
    
    arr = fread(fid,inf,'double');
    n = arr(1);
    x = arr(2:n+1)';
    F = arr(n+2:2*n+1)';    
    fclose(fid);     
end