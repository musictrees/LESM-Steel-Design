%% LESM launcher function
%
% This function is called to initialize the LESM program.
% It provides two options for running the program: a graphical and a
% nongraphical mode.
%
function launch(mode,fileName)
    if nargin == 0
        mode = 1;
    end
    
    warning('off');
    
    % Close current instance
    if (strcmp(class(findall(0)),'matlab.graphics.Graphics')) %#ok<STISA>
        delete(findall(0));
    end
    
    % Switch between nongraphical and graphical modes
    if mode == 0 % NONGRAPHICAL MODE
        % Initialize model object
        model = Model();
        
        % Open input file
        fid = fopen(fileName,'rt');
        
        % Check for valid input file
        if fid > 0
            % Read model information
            fprintf(1,'Pre-processing...\n');
            [vs,print,~,lc] = readFile(fid,model,mode);
            
            % Check input file version compatibility
            if vs == 1
                % Check analysis option and create analysis driver object
                include_constants;
                switch model.whichSolver
                    case STATIC_LINEAR
                        model.drv = Drv_LES(mode,model);
                    case DYNAMIC_NEWMARK_LINEAR
                        model.drv = Drv_LED(DYNAMIC_NEWMARK_LINEAR,mode,model);
                    case DYNAMIC_MODALSUP_LINEAR
                        model.drv = Drv_LED(DYNAMIC_MODALSUP_LINEAR,mode,model);
                end
                
                % Process provided data
                status = model.drv.process();
                
                % Print analysis results
                if status == 1
                    if lc <= model.nlc
                        string = model.strLc{lc};
                    else
                        string = model.strComb{lc-model.nlc};
                    end
                    print.results(lc,string);
                end
            else
                fprintf(1,'This file version is not compatible with the program version!\n');
            end
        else
            fprintf(1,'Error opening input file!\n');
        end
        
        % Clear memory
        clear;
        
    elseif mode == 1 % GRAPHICAL MODE
        % Show splashscreen
        %if ~isdeployed
        %    s = SplashScreen('LESM','logo_lesm_splash.png');
        %end
        
        % Build main graphical interface dialog
        GUI_Main;
        
        % Delete splashScreen
        %if ~isdeployed
        %    delete(s);
        %end
        
        % Make GUI visible
        gui = findall(groot,'Type','figure');
        set(gui,'Visible','on');
        
    else % INVALID MODE
        fprintf(1,'Inavlid mode!\n');
        fprintf(1,'mode = 0 -> Nongraphical\n');
        fprintf(1,'mode = 1 -> Graphical\n');
        clear;
    end
end