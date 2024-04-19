%% Read file Function
%
% This is an auxiliary file that contains functions to read a
% neutral-format file with the _.lsm_ extension 
% This neutral-format file contains all information about a linear elements
% structural model.
%
%% Main function
% Output:
%  vs:       flag for version compatibility
%  print:    object of the Print class
%  draw:     object of the Draw class
%  nclc:     current load case id
% Input arguments:
%  fid:      integer identifier of the input file
%  model:    handle to an object of the Model class
%  GUI_Mode: flag for the graphical version of LESM being used
function [vs,print,draw,nclc] = readFile(fid,model,GUI_Mode,pathname)
    if nargin == 2
        GUI_Mode = true;
        pathname = '';
    elseif nargin == 3
        pathname = '';
    end
    
    % Reads input file
    next_line = true;
    while next_line && ~feof(fid)
        % Get file line
        tline = fgetl(fid);
        % Get rid of blank spaces
        string = deblank(tline);
        
        % Look for for tag strings
        switch string
            case '%HEADER.VERSION'
                [vs,version] = readVersion(fid);
                if (version == 1.0) || (version == 1.1)
                    [print,draw,nclc] = readFile_v1(fid,model,GUI_Mode); 
                elseif version == 2.0
                    [print,draw,nclc] = readFile_v2(fid,model,GUI_Mode);
                elseif version == 2.1
                    [print,draw,nclc] = readFile_v2_1(fid,model,GUI_Mode);
                elseif version == 3.0
                    [print,draw,nclc] = readFile_v3(fid,model,GUI_Mode,pathname);
                else
                    print = [];
                    draw = [];
                    nclc = 0;
                end
                next_line = false;

            case '%END'
                vs = 0;
                print = [];
                draw = [];
                nclc = 0;
                next_line = false;
        end
    end
    
    % Close file
    fclose(fid);
end

%% Auxiliary functions
%--------------------------------------------------------------------------
function [vs,version] = readVersion(fid)
    % Read program version number
    version = fscanf(fid,'%f',1);
    
    % Check version compatibility
    if (version == 1.0) || (version == 1.1)  || (version == 2.0) ||...
       (version == 2.1) || (version == 3.0)
        vs = 1;
    else
        vs = 0;
        return
    end
end