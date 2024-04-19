classdef NBR8800<  handle & Standards
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
   
    
   
   properties
       ActionList={'Permanente';'Variável'};
        Type1List={'Normais';'Especiais ou de contrução';'Excepcionais'};
        Type2PermanentList={'Peso própio de estr.metálica';'Peso própio de estr.pré-moldada';'Peso própio de estr.moldadas no local/elem.construtivos industrializados e empuxos permanentes';'Peso própio de elem.contrutivos industrializados com adições in loco';...
            'Peso própio de elem.construtivos  em geral e equipamentos';'Indiretas'};
        Type2VariableList={'Efeito Temperatura';'Ação Vento';'Acões Truncadas';'Outros'};
        
        Type2PsiList={'Locais em que não há predominancia de pesos e de equipamentos que permanecem fixos por longos periodos,nem de elevadas concentrações de pessoas';
            'Locais em que há predominancia de pesos e de equipamentos que permanecem fixos por longos periodos,ou de elevadas concentrações de pessoas';
            'Bibliotecas,arquivos, depósitos,oficinas,garagens e sobrecargas de temperatura';
            'Pressão dinâmica do vento nas estruturas em geral';
            'variações uniformes de temperatura em relação á média anual local';
            'Passarelas de pedestres';
            'Vigas de rolamento de pontes rolantes';
            'Pilares e outros elementos ou subestruturas que suportam vigas de rolamento de pontes rolantes';%REFAZER
            }
        
        
        permanentValues=[1.25;1.30;1.35;1.40;1.50;1.20;1.15;1.20;1.25;1.30;1.40;1.20;1.10;1.15;1.15;1.20;1.30;0];
        keyPermanentValues={'1,1,1';'1,1,2';'1,1,3';'1,1,4';'1,1,5';'1,1,6';'1,2,1';'1,2,2';'1,2,3';'1,2,4';'1,2,5';'1,2,6';'1,3,1';'1,3,2';'1,3,3';'1,3,4';'1,3,5';'1,3,6'}
        
        variableValues=[1.20;1.40;1.20;1.50;1.0;1.20;1.10;1.30;1.0;1.0;1.0;1.0];
        keyVariableValues={'2,1,1';'2,1,2';'2,1,3';'2,1,4';'2,2,1';'2,2,2';'2,2,3';'2,2,4';'2,3,1';'2,3,2';'2,3,3';'2,3,4'}
        
        PsiValues=[0.5;0.7;0.8;0.6;0.6;0.6;1.0;0.7];%REFAZER
        keyPsiValues={'1';'2';'3';'4';'5';'6';'7';'8'}%REFAZER
        
        
        
        hashTablePermanent=0;
        hashTableVariable=0;
        hashTablePsi=0;
        designReportFile;
        
        
   end
   
    
    
    methods
        
        function Standards=NBR8800(reportFile)
            Standards.designReport=DesignReport();
                         Standards.hashTablePermanent=containers.Map(Standards.keyPermanentValues,Standards.permanentValues);
                         Standards.hashTableVariable=containers.Map(Standards.keyVariableValues,Standards.variableValues);
                         Standards.hashTablePsi=containers.Map(Standards.keyPsiValues,Standards.PsiValues);
            Standards.designReportFile=reportFile;
            
            
      
        end
        
        function coeficient = LoadCombination(Standards,keyValue,psi)
            if nargin==2
                x=split(keyValue,",");
                if x(1)=="1"
                    coeficient = Standards.hashTablePermanent(keyValue);
                else
                    coeficient = Standards.hashTableVariable(keyValue);
                end
            else
                coeficient=Standards.hashTablePsi(keyValue);
            end
            
        end
        
        %Principal Function
        function solve(Standards,element,designReportFile,i)
          
         
          
            %             element.section.mpl
            %             element.section.Wx*element.material.leakage*0.7
           
            %             Standards.fltSolver(element,fid,i);
            %Standards.flmSolver(element);
            %Standards.flaSolver (element);
             section=[  "<div id='Compres-" + i + "' class='card border-titledark col-md-12 bg-light mt-5' style='padding: 0px;'>" 
             "<div class='card-header font-weight-bold'>elemento</div>" 
              
              "<div class='row align-items-center  p-3'>" 
               "<div class='col-md-12'>" 
               "<h5 class='card-title'>Elemento-" + i + "</h5>" 
               "</div>"
               "<div class='card-body '>" %%subistituir por função
               "<div class='row'>"
                "<div class='col-md-6 text-right'>" 
                   "<p class='font-weight-bold'>"+"Material Infos</p>"
            
                   "<p>"+"<span class='font-weight-bold'>E</span>="+element.material.elasticity+"</p>"
                   "<p>"+"<span class='font-weight-bold'>G</span>="+element.material.shear+"</p>"
                   "<p>"+"<span class='font-weight-bold'>Fy</span>="+element.material.leakage+"</p>"
                   "<p>"+"<span class='font-weight-bold'>Fu</span>="+element.material.yield+"</p>"
                   
                "</div >" 
                "<div class='col-md-6 text-right'>" 
                    "<img src='img_ieh_laminados.gif'/>"
                "</div >" 
              "</div >" 
              "</div >" 
             
              "</div></div>" %close card
                
             ];
            fprintf(designReportFile,'%s',section);
                 %compressao
                %fprintf(designReportFile,'%s',Standards.designReport.openSection(i,"Compress"));
              
% %             disp(Standards.getNexFactor(element,Standards.getKfactor(element)));
% %              disp(Standards.getNeyFactor(element,Standards.getKfactor(element)));
% %               disp("QA"+ Standards.getLocalBucklingFactorQa(element,4));
% %                disp( "QS"+Standards.getLocalBucklingFactorQs(element,4));
%                Q=(Standards.getLocalBucklingFactorQa(element,4)*Standards.getLocalBucklingFactorQs(element,4));


%               X= Standards.getQsiFactor(element,Q,Ne);
%              
%              Nrd=(X*Q*element.section.Ag*element.material.leakage)/(1.1)
           
            %flexão
            fprintf(designReportFile,'%s', Standards.designReport.openCardSection(i,"traçao",element.maxTraction,element.maxTractionId));
            Standards.NrdSolve(element,designReportFile,i); %TRAÇAO FUNÇAO
            closeSection=["</div></div></div>"];
            fprintf(designReportFile,'%s',closeSection);
          
             fprintf(designReportFile,'%s', Standards.designReport.openCardSection(i,"flex-x-x",element.maxBendMoment_XY,element.maxBendMoment_XY_CombinationId));
            
                        [classification,lambda,lambdaP,lambdaR,Mpl,Mr,Mcr]=Standards.classifyFlm(element,1,designReportFile);
                        
                        MrdFLM=Standards.MrdFlmFlaProcess(classification,lambda,lambdaP,lambdaR,Mpl,Mr,Mcr,"FLM",designReportFile,element,i);
                        [classification,lambda,lambdaP,lambdaR,Mpl,Mr,Mcr]=Standards.classifyFla(element,1,designReportFile);
                        MrdFLA=Standards.MrdFlmFlaProcess(classification,lambda,lambdaP,lambdaR,Mpl,Mr,Mcr,"FLA",designReportFile,element,i);
                        [classification,lambda,lambdaP,lambdaR,Mpl,Mr,Mcr,cb]=Standards.classifyFlt(element,1,designReportFile);
                        MrdFLT=Standards.MrdFltProcess(classification,lambda,lambdaP,lambdaR,Mpl,Mr,Mcr,cb,designReportFile,element,i);
%             
                closeSection=["</div></div></div>"];
              fprintf(designReportFile,'%s',closeSection);
%             Standards.getKfactor(element)
            
            
%            
               %compress

            fprintf(designReportFile,'%s',Standards.designReport.openCardSection(i,'Compress',element.maxCompression,element.maxCompressionId));
            [Kx,Ky]=Standards.getKfactor(element);
           
            Ne=Standards.getNeValue(Standards.getNexFactor(element,Kx),Standards.getNeyFactor(element,Ky),Standards.getNezFactor(element));
            Qa=Standards.getLocalBucklingFactorQa(Standards,element,4);
            Qs=Standards.getLocalBucklingFactorQs(Standards,element,4);
            Q=(Qa*Qs);
            fprintf(designReportFile,'%s',Standards.designReport.QSection(Qa,Qs));
            X=Standards.getQsiFactor(Standards,element,Q,Ne);
            Nrd=CompressSolve(Standards,X,Q,element.section.Ag,element.material.leakage,i);
            fprintf(designReportFile,'%s', "</div></div></div>");  %close card
            
              %shear
            fprintf(designReportFile,'%s',Standards.designReport.openCardSection(i,'Shear',element.maxShearForce_XY,element.maxShearForce_XY_CombinationId));
            Standards.shearResistence(element,1,designReportFile,i);
            fprintf(designReportFile,'%s', "</div></div>");  %close card
            
            Standards.axialBendingCombination(element);
        end
        %FLEXAO SIMPLES
        function [classification,lambda,lambdaP,lambdaR,Mpl,Mr,Mcr]=classifyFlm(Standards,element,geometricType,designReportFile)
            switch geometricType
                case 1
                    disp('MRY');
                    mry=(element.material.leakage-element.material.leakage*0.3)*element.section.Wy
                    mrySemRed=(element.material.leakage)*element.section.Wy
                    Mpl=element.section.Zx*element.material.leakage;%REVISAR CHECAR VALOR
                    lambda=(element.section.bf/2)/element.section.tf
                    lambdaP=0.38*sqrt(element.material.elasticity/element.material.leakage)
                    lambdaR=0.83*sqrt(element.material.elasticity/(element.material.leakage-element.material.leakage*0.3))
                    Wc=(element.section.Ix/element.section.d/2); %%Onde dy é a distância do centro de gravidade até a borda da mesa comprimida para esse caso de viga i esata dy=d/2 de maneira provisoria
                    Mr=(element.material.leakage-element.material.leakage*0.3)*element.section.Wx     %%REVISAR WX saber em qual momento usar wx ou wy
                    Mcr=((0.69*element.material.elasticity)/pi^2)*Wc; %%Revisar  Mcr pagina 136 topico 6 laminado Wc=ix/dy ONDE :Onde dy é a distância do centro de gravidade até a borda da mesa comprimida para esse caso de viga i esata dy=d/2 de maneira provisoria
                     data.Lambda="$$\lambda=\frac{b}{t}$$";
                     data.LambdaP="$$\lambda_p=0.38*\sqrt{\frac{E}{f_y}}$$";
                     data.LambdaR="$$\lambda_r=0.83*\sqrt{\frac{E}{f_y-\sigma_r}}$$";
                     data.title="Verificação Flambagem local da Mesa (FLM)";
                     fprintf(designReportFile,'%s',Standards.designReport.flexSection(lambda,lambdaP,lambdaR,Mr,Mcr,1,data));
                    
                case 2
                    Mpl=element.section.mpl;
                    lambda=(element.section.bf/2)/element.section.tf;
                    lambdaP=0.38*sqrt(element.material.elasticity/element.material.leakage);
                    lambdaR=0.83*sqrt(element.material.elasticity/(element.material.leakage-element.material.leakage*0.3));
                    Wc=(element.section.Ix/element.section.d/2); %%Onde dy é a distância do centro de gravidade até a borda da mesa comprimida para esse caso de viga i esata dy=d/2 de maneira provisoria
                    Mr=(element.material.leakage-element.material.leakage*0.3)*Wc;
                    Mcr=((0.69*element.material.elasticity)/pi^2)*Wc;
                    
                case 3
                    Wc=(element.section.Ix/element.section.d/2); %%Onde dy é a distância do centro de gravidade até a borda da mesa comprimida para esse caso de viga i esata dy=d/2 de maneira provisoria SUBSTITUIR IX POR IY?
                    Mpl=element.section.mply;%conferir formula
                    Mr=(element.material.leakage-element.material.leakage*0.3)*element.section.Wy;%conferir se usa Wy
                    Mcr=((0.69*element.material.elasticity)/pi^2)*Wc;
                    
                    lambda=(element.section.bf/2)/element.section.tf;  %tirar duvida se as medidas sao iguais para eixo fraco %REVISAR
                    lambdaP=0.38*sqrt(element.material.elasticity/element.material.leakage);
                    lambdaR=0.83*sqrt(element.material.elasticity/(element.material.leakage-element.material.leakage*0.3));
                    
                case 4
                    result="nao possui FLM"
                    
                case 5
                    %completar
                    
                otherwise
            end
            if(lambda <= lambdaP)
                classification="compacta";
                return
                
            end
            
            if ((lambdaP<lambda) && (lambda<=lambdaR))
                classification="semiCompacta";
                return
            end
            
            if (lambda>lambdaR)
                classification="esbelta" ;
                return
            end
            
            
        end
        
        function [classification,lambda,lambdaP,lambdaR,Mpl,Mr,Mcr]=classifyFla(Standards,element,geometricType,designReportFile) %Tipo de seção e eixo de flexão pagina 134    TABELA G.1 PARA FLM
            
            switch geometricType
                case 1
                    disp('valor de Zy');
                    disp(element.section.Zy);
                    Mpl=element.section.Zx*element.material.leakage;
                    Mr=element.material.leakage*element.section.Wx;
                    Mcr=0;%ANEXO H REVISAR COM ANTONIONE
                    lambda=element.section.h/element.section.tw;
                    lambdaP=3.76*sqrt((element.material.elasticity/element.material.leakage));
                    lambdaR=5.70*sqrt((element.material.elasticity/element.material.leakage));
                     data.Lambda="$$\lambda=\frac{h}{tw}$$";
                     data.LambdaP="$$\lambda_p=3.76*\sqrt{\frac{E}{f_y}}$$";
                     data.LambdaR="$$\lambda_r=5.70*\sqrt{\frac{E}{f_y-\sigma_r}}$$";
                     data.title="Verificação Flambagem local da Alma (FLA)";
                     fprintf(designReportFile,'%s',Standards.designReport.flexSection(lambda,lambdaP,lambdaR,Mr,Mcr,1,data));
                case 2
                    Mpl=element.section.Zx*element.material.leakage;
                    Hc;%calcular Hc
                    Hp;
                    
                    Mr=element.material.leakage*element.section.Wx;
                    Mcr;%ANEXO H REVISAR COM ANTONIONE
                    lambda;%hc/tw
                    lambdaP=((Hc/Hp)*sqrt((element.material.elasticity/element.material.leakage)))/((0.54*(Mpl/Mr)-0.09))^2;
                    lambdaR=5.70*sqrt((element.material.elasticity/element.material.leakage));
                    
                case 3
                    
                    Mpl=element.section.Zy*element.material.leakage;
                    Mr;%completar
                    Mcr;%completar
                    lambda;%
                    lambdaP=1.12*sqrt((element.material.elasticity/element.material.leakage));
                    lambdaR=1.40*sqrt((element.material.elasticity/element.material.leakage));
                case 4
                    result="nao possui FLA"
                    
                case 5
                    %completar
                    
                otherwise
            end
            if(lambda <= lambdaP)
                classification="compacta";
                return
                
            end
            
            if ((lambdaP<lambda) && (lambda<=lambdaR))
                classification="semiCompacta";
                return
            end
            
            if (lambda>lambdaR)
                classification="esbelta" ;
                return
            end
            
            
        end
        
        
         %funçao aux para achar cb essa funcao retorna a posiçao local do
         %elemento que desejo o momento partindo da divisao global e a do
         %numero do elemento correspondente.
        
        function normalizeInfos=normalizeLength(~,element,division)
            piece=element.length/division;
            
            sumOutLength=0;
            for i=1:length(element.elems)
                element.elems(i).length;
                if((sumOutLength+element.elems(i).length)<piece)
                    sumOutLength=sumOutLength+element.elems(i).length;
                    
                else
                    elementIdNumber=i;
                    momentLocalPosition=piece-sumOutLength;
                    normalizeInfos=[elementIdNumber,momentLocalPosition];
                    break
                end
            end
        end
        function [classification,lambda,lambdaP,lambdaR,Mpl,Mr,Mcr,cb]=classifyFlt(Standards,element,geometricType,designReportFile) %Tipo de seção e eixo de flexão pagina 134    TABELA G.1 PARA FLM
            
            switch geometricType
                case 1
%                     [M,~] = element.intBendingMoment_XY([(element.length/4),(element.length/2),((element.length*3)/4)]);

            position1=Standards.normalizeLength(element,4);
            position2=Standards.normalizeLength(element,2);
            position3=Standards.normalizeLength(element,(4/3));
        
            [M1,~] = element.elems(position1(1)).intBendingMoment_XY((position1(2)));
             [M2,~] = element.elems(position2(1)).intBendingMoment_XY((position2(2)));
              [M3,~] = element.elems(position3(1)).intBendingMoment_XY((position3(2)));

                M=cat(2,M1,M2,M3)

                    Mpl=element.section.Zx*element.material.leakage;
                    disp('CB AQUI');
%                     cb = element.maxBendMoment_XY(1) *12.5 /(element.maxBendMoment_XY(1)*2.5 + M(1)*3 + M(2)*4 + M(3)*3)%PERGUNTAR PARA O PEDRO A EXTREMIDADE QUE COMEÇA A DIVISAO e com pegar exato pag 47 nbr
                    cb=element.Cb %MUDADO LENGTH PARA Cb
                    disp((element.section.Wx));
                    disp(element.material.leakage);
                    disp(element.section.J);
                    b1 = (0.7 * element.section.Wx*element.material.leakage)/((element.material.elasticity * element.section.J)) %PAG135  element.section.J subs por 1.72x10-8CONFERIR FORMULA ADE J e B1
                    Mr=element.material.leakage*element.section.Wx;
                    Mcr=((cb * pi^2 * element.material.elasticity * element.section.Iy) / element.length^2) * sqrt(element.section.cw/element.section.Iy * (1 + 0.039 * element.section.J * element.length^2/element.section.cw) );%element.section.J subs por 1.72x10-8CONFERIR FORMULA ADE J e para Lb estouy usando Lenght provisoriamente
                    lambda=element.Lb/element.section.ry; %MUDADO LENGTH PARA LB
                    lambdaP=1.76*sqrt(element.material.elasticity/element.material.leakage);
                    lambdaR=1.38*sqrt(element.section.Iy * element.section.J)/(element.section.ry * element.section.J * b1) * sqrt( 1 + sqrt( 1 + (27*element.section.cw * b1^2 / element.section.Iy)));%element.section.J subs por 1.72x10-8CONFERIR FORMULA ADE J
                     data.Lambda="$$\lambda=\frac{Lb}{r_y}$$";
                     data.LambdaP="$$\lambda_p=1.76*\sqrt{\frac{E}{f_y}}$$";
                     data.LambdaR="$$\lambda_r=\frac{1.38*\sqrt{I_y*J}}{r_y*J*\beta _1}*\sqrt{1+\sqrt{}\frac{27*Cw*\beta_1^{2}}{I_y}}$$";
                     data.title="Verificação Flambagem lateral a torção (FLT)";
                     fprintf(designReportFile,'%s',Standards.designReport.flexSection(lambda,lambdaP,lambdaR,Mr,Mcr,1,data));
                  
                case 2
                    
                    
                case 3
                    
                case 4
                    
                    
                case 5
                    
                    
                otherwise
            end
            if(lambda <= lambdaP)
                classification="compacta";
                return
                
            end
            
            if ((lambdaP<lambda) && (lambda<=lambdaR))
                classification="semiCompacta";
                return
            end
            
            if (lambda>lambdaR)
                classification="esbelta" ;
                return
            end
            
            
        end
        
        function MrdFlmFla=MrdFlmFlaProcess(Standards,classification,lambda,lambdaP,lambdaR,Mpl,Mr,Mcr,flagInfo,designReportFile,element,i)
           
            switch classification
                case "compacta"
                    MrdFlmFla=Mpl/1.1;
                    data.MrdString="$$Mrd=\frac{M_pl}{\gamma_a1}$$";
                    data.Mrd=MrdFlmFla;
                case "semiCompacta"
                    MrdFlmFla=(1/1.1)*(Mpl-(Mpl-Mr)*((lambda-lambdaP)/(lambdaR-lambdaP)));%substituir 1.1 por uma propreidade da norma
                   
                    data.MrdString="$$Mrd=\frac{1}{\gamma_a1}*[M_p1-(M_p1-M_r)*\frac{\lambda-\lambda_p}{\lambda_r - \lambda_p})]$$"
                    data.Mrd=MrdFlmFla;
                case "esbelta"
                    if(flagInfo=="FLA")
                        %completar com anexo H
                    else
                        MrdFlmFla=Mcr/1.1
                    end
                    
                otherwise
            end
           data.Msd=element.maxBendMoment_XY;
            if(MrdFlmFla>=(element.maxBendMoment_XY*1000))
               
                data.approves=true;
            else
                data.approves=false;
            end
            allStructuralGroups=getappdata(0,'allStructuralGroups');
           %%
            if length(allStructuralGroups(i).bendingXYVerify)==3%%apagar para voltar antes
                %%%%apagar
                if(flagInfo=="FLA")
                    verifyInfo.type='FLA';
                    verifyInfo.rd=MrdFlmFla;
                    allStructuralGroups(i).bendingXYVerify(2)=verifyInfo;
                else
                    verifyInfo.type='FLM';
                    verifyInfo.rd=MrdFlmFla;
                    allStructuralGroups(i).bendingXYVerify(1)=verifyInfo;
                end
            else %%apagar
                if(flagInfo=="FLA")
                    verifyInfo.type='FLA';
                    verifyInfo.rd=MrdFlmFla;
                    allStructuralGroups(i).bendingXYVerify=[allStructuralGroups(i).bendingXYVerify,verifyInfo];
                else
                    verifyInfo.type='FLM';
                    verifyInfo.rd=MrdFlmFla;
                    allStructuralGroups(i).bendingXYVerify=[allStructuralGroups(i).bendingXYVerify,verifyInfo];
                end
            end%%apagar
           
          
            %UPDATE STRUCTURALGROUPS AND MODEL
            setappdata(0,'allStructuralGroups',allStructuralGroups);
            model = getappdata(0,'model');
            model.allStructuralGroups=allStructuralGroups;
            setappdata(0,'model',model);
            fprintf(designReportFile,'%s', Standards.designReport.flexMrd(data));
            %UPDATE STRUCTURALGROUPS AND MODEL
            
            
        end
        
        function MrdFlt=MrdFltProcess(Standards,classification,lambda,lambdaP,lambdaR,Mpl,Mr,Mcr,cb,designReportFile,element,i)
           disp('TESTE');
             Mcr
            switch classification
                case "compacta"
                     MrdFlt=Mpl/1.1;
                     data.MrdString="$$Mrd=\frac{M_pl}{\gamma_a1}$$";
                     data.Mrd=MrdFlt;
                   
                case "semiCompacta"
                    MrdFlt=(cb/1.1) * (Mpl - (Mpl -  Mr)*((lambda - lambdaP)/(lambdaR - lambdaP)));
                    data.MrdString="$$Mrd=\frac{C_b1}{\gamma_a1}*[M_p1-(M_p1-M_r)*\frac{\lambda-\lambda_p}{\lambda_r - \lambda_p})]$$";
                    data.Mrd=MrdFlt;
                   
                case "esbelta"
                    MrdFlt=Mcr/1.1;
                     data.MrdString="$$Mrd=\frac{M_cr}{\gamma_a1}$$";
                     data.Mrd=MrdFlt;
                otherwise
                     
            end
               allStructuralGroups=getappdata(0,'allStructuralGroups');
            
            if length(allStructuralGroups(i).bendingXYVerify)==3%%apagar para voltar antes
                %%%%apagar
               
                    verifyInfo.type='FLT';
                    verifyInfo.rd=MrdFlt;
                    allStructuralGroups(i).bendingXYVerify(3)=verifyInfo;
             
            else %%apagar
                    verifyInfo.type='FLT';
                    verifyInfo.rd=MrdFlt;
                     allStructuralGroups(i).bendingXYVerify=[allStructuralGroups(i).bendingXYVerify,verifyInfo];
             
            end%%apagar
           
              if(MrdFlt>=element.maxBendMoment_XY*1000)
                data.approves=true;
                disp("aqui")
              else
                disp("aqui2")
                data.approves=false;
            end
              fprintf(designReportFile,'%s', Standards.designReport.flexMrd(data));

             
        end
        
        
        
        
        %%COMPRESSAO
        %funçao para pegar o valor de k
        function [Kx,Ky]= getKfactor(Standards,element)
%             fixed=[1,1,0,0,0,1]; %node class olhar
%             pinned=[1,1,0,0,0,0];
%             rotFixedTransFree=[0,0,0,0,0,1];
%             free=[0,0,0,0,0,0];
%             data.boundary="";
%            
%             if(element.nodes(1).ebc==fixed & element.nodes(2).ebc==fixed)
%                 k=0.5;
%                 data.boundary="fixed";
%             end
%             if ((element.nodes(1).ebc==pinned & element.nodes(2).ebc==pinned))
%                  k=1;
%                  data.boundary="pinned";
%             end
%             if ((element.nodes(1).ebc==fixed & element.nodes(2).ebc==free) | (element.nodes(1).ebc==free & element.nodes(2).ebc==fixed))
%                  k=2;
%                    data.boundary="free/fixed";
%             end
%             
%             if ((element.nodes(1).ebc==pinned & element.nodes(2).ebc==rotFixedTransFree) | (element.nodes(1).ebc==rotFixedTransFree & element.nodes(2).ebc==pinned))
%                  k=2;
%                   data.boundary="pinned/rotFixedTransFree";
%             end
%             
%             if ((element.nodes(1).ebc==fixed & element.nodes(2).ebc==pinned) | (element.nodes(1).ebc==pinned & element.nodes(2).ebc==fixed))
%                  k=0.7;
%                  data.boundary="fixed/pinned";
%             end
            Kx=element.Kx;%valor mudado
            Ky=element.Ky;%valor mudado
            data.KxValue=Kx;
            data.KyValue=Ky;
            
             fprintf(Standards.designReportFile,'%s',Standards.designReport.kFactorSection(data));
            
        end
        
        function Nex= getNexFactor(~,element,Kx)
            Nex.value=((pi^2)*element.material.elasticity*element.section.Ix)/(Kx*element.Lx)^2;
            Nex.equationString="$$\frac{\pi^{2}*E*I_{x}}{(k_{x}*L{_{x})}^{2}}$$";
            Nex.equationValue="$$\frac{\pi^{2}*"+element.material.elasticity+"*"+element.section.Ix+"}{("+Kx+"*"+element.length+")}^{2}$$";
            
        end
        function Ney=getNeyFactor(~,element,Ky)
            Ney.value=((pi^2)*element.material.elasticity*element.section.Iy)/(Ky*element.Ly)^2;
            
            Ney.equationString="$$\frac{\pi^{2}*E*I_{y}}{(k_{x}*L{_{y})}^{2}}$$";
            Ney.equationValue="$$\frac{\pi^{2}*"+element.material.elasticity+"*"+element.section.Iy+"}{("+Ky+"*"+element.length+")}^{2}$$";
        end
        function Nez= getNezFactor(~,element)
         
            Nez.equationString="$$\frac{1}{r_0^{2}}[\frac{\pi^{2}ECw}{(KzLz)^{2}}+GJ]$$";
            Nez.equationValue="$$\frac{1}{"+element.section.r0+"^{2}}[\frac{\pi^{2}"+element.material.elasticity+"*"+element.section.cw+"}{("+element.material.shear+""+element.Lb+")^{2}}+"+element.material.shear+"*"+element.section.J+"]$$";
            Nez.value=(1/(element.section.r0^2))*((((pi^2)*element.material.elasticity*element.section.cw)/((element.Kz*element.Lb)^2))+(element.material.shear*element.section.J));%element.Kz*element.Lb mudado
        end
        
        function Ne= getNeValue(Standards,Nex,Ney,Nez)
            
            fprintf(Standards.designReportFile,'%s',Standards.designReport.NeSection(Nex,Ney,Nez));
            %put Nz in min function !IMPORTANT
            Ne=min([Nex.value,Ney.value,Nez.value]);
        end
        
        
        
        function  Qa=getLocalBucklingFactorQa(~,Standards,element,geometricType)
            
             
            switch geometricType
                case 4
                    esbeltezLim.value=1.49*(sqrt((element.material.elasticity/element.material.leakage))); %tabela g.1 pag 128 BOTAR PRA FORA E TESTAR
                    esbeltez.value=(element.section.bw)/element.section.tw;
                    esbeltez.equationString="$$\frac{b_{w}}{t_{w}}$$";
                    esbeltez.equationValue="$$\frac{"+element.section.bw+"}{"+element.section.tw+"}$$";
                    esbeltezLim.equationString="$$1.49*\sqrt{\frac{E}{F_{y}}}$$";
                    esbeltezLim.equationValue="$$1.49*\sqrt{\frac{"+element.material.elasticity+"}{"+element.material.leakage+"}}$$";
                    specificCase.befEquationString=0;
                    
                    if(esbeltez.value <= esbeltezLim.value )
                        Qa=1;
                        specificCase.qaValue="$$Qa="+Qa+"$$";
                     
                    else
                        
                        bef=1.92*element.section.tw*sqrt((element.material.elasticity/element.material.leakage))*(1-(0.34/(element.section.bw/element.section.tw)*(sqrt((element.material.elasticity/element.material.leakage)))));
                       
                        Aef=element.section.Ag-((element.section.bw-bef)*element.section.tw);
                        Qa=(Aef/element.section.Ag);
                        
                        
                        
                        specificCase.befEquationString="$$Bef=1.92*t_{w}*\sqrt{\frac{E}{f_{y}}}*[1-\frac{C_{a}}{\frac{b}{t}}*\sqrt{\frac{E}{f_{y}}}]$$";
                        specificCase.befEquationValue="$$Bef=1.92*"+element.section.tw+"*\sqrt{\frac{"+element.material.elasticity+"}{"+element.material.leakage+"}}*[1-\frac{"+0.38+"}{\frac{"+element.section.bw+"}{"+element.section.tw+"}}*\sqrt{\frac{"+element.material.elasticity+"}{{"+element.material.leakage+"}}}]$$";
                        specificCase.befValue="$$Bef="+bef+"$$";
                        
                        specificCase.aefEquationString="$$Ag=A_{g}-\sum_{}^{}((b_{w}-b_{ef})*t_{w})$$";
                        specificCase.aefValue="$$Aef="+Aef+"$$";
                        
                        specificCase.qaEquationString="$$Qa=\frac{A_{ef}}{A_{g}}$$";
                        specificCase.qaEquationValue="$$Qa=\frac{"+Aef+"}{"+element.section.Ag+"}$$";
                        specificCase.qaValue="$$Qa="+Qa+"$$";
                    end
                 
                    
                    
                  
                    
                    fprintf(Standards.designReportFile,'%s',Standards.designReport.QaSection(esbeltez,esbeltezLim,specificCase));
                    
                    
                otherwise
            end
           
           
        end
        
        function Qs= getLocalBucklingFactorQs(~,Standards,element,geometricType)
            esbeltez.value=(element.section.bf/2)/element.section.tf;
            esbeltez.equationString="$$\frac{\frac{b_{f}}{2}}{tf}$$";
            esbeltez.equationValue="$$\frac{\frac{"+element.section.bf+"}{2}}{"+element.section.tf+"}$$";
            
            specificCase.qsValue=0;
            switch geometricType
             case 4
                   esbeltezLim.value=0.56*sqrt((element.material.elasticity/element.material.leakage));
                   esbeltezLim.equationString="$$0.56*\sqrt{\frac{E}{f_{y}}}$$";
                   if(mean([esbeltez.value]) <= mean([esbeltezLim.value]))
                       Qs=1;
                       specificCase.qsValue="$$Qs="+Qs+"$$";
                       specificCase.compare="$$"+esbeltez.value+"<="+esbeltezLim.value+"$$";
                       specificCase.compareCase=1;
                   
                           
                   else
                       if(esbeltezLim <  esbeltez.value && esbeltez.value < (1.03*sqrt((element.material.elasticity/element.material.leakage))))
                           Qs=1.415-0.74*(esbeltez)*sqrt((element.material.leakage/element.material.elasticity));
                            specificCase.qsValue="$$Qs="+Qs+"$$";
                            specificCase.qsEquationString="$$Qs=1.415-0.74*\frac{b}{t}*\sqrt{\frac{f_{y}}{E}}$$";
                            specificCase.qsEquationValue="$$Qs=1.415-0.74*"+esbeltez.value+"*\sqrt{\frac{"+element.material.leakage+"}{"+element.material.elasticity+"}}$$";
                            specificCase.esbeltezLimSup="$$Qs=1.03\sqrt{\frac{E}{f_{y}}}$$";
                            specificCase.esbeltezLimSupValue=(1.03*sqrt((element.material.elasticity/element.material.leakage)));
                            specificCase.compare="$$"+esbeltezLim.value+"<"+ esbeltez.value+"<=" + specificCase.esbeltezLimSupValue+ "$$";
                          specificCase.compareCase=2;
                       end
                       if(  esbeltez.value > 1.03*sqrt((element.material.elasticity/element.material.leakage)))
                              
                       Qs=0.69*element.material.elasticity/element.material.leakage*((esbeltez.value)^2);
                           specificCase.qsValue="$$Qs="+Qs+"$$";
                           specificCase.qsEquationString="$$Qs=\frac{0.69*E}{f_{y}*(\frac{b}{t})^{2}}$$";
                           specificCase.qsEquationValue="$$Qs=\frac{0.69*"+element.material.elasticity+"}{"+element.material.elasticity+"*("+esbeltez.value+")^{2}}$$";
                           specificCase.esbeltezLimSup="$$Qs=1.03\sqrt{\frac{E}{f_{y}}}$$";
                           specificCase.esbeltezLimSupValue=(1.03*sqrt((element.material.elasticity/element.material.leakage)));
                           specificCase.compare="$$"+esbeltez.value+">" + specificCase.esbeltezLimSupValue+ "$$";
                           specificCase.compareCase=3;
                       end
                   end
             otherwise
                    
            end
            
            
            
            
  fprintf(Standards.designReportFile,'%s',Standards.designReport.QsSection(esbeltez,esbeltezLim,specificCase));

            
        end
        
        function X = getQsiFactor(~,Standards,element,Q,Ne)
           
            esbeltezReduzida.equationString="$$\gamma_{o}=\sqrt{Q*Ag*\frac{f_{y}}{Ne}}$$";
             esbeltezReduzida.equationValue="$$\gamma_{o}=\sqrt{"+Q+"*"+element.section.Ag+"*\frac{"+element.material.leakage+"}{"+Ne+"}}$$";
              esbeltezReduzida.value=sqrt(Q*element.section.Ag*element.material.leakage/Ne);
              
            if(esbeltezReduzida.value<= 1.5)
                X=0.658^(esbeltezReduzida.value^2);
                specificCase.compareCase=1;
                specificCase.compare="$$\gamma_{0}<=1.5$$"
                specificCase.equationString="$$\chi=0.648^{\gamma_{0}^{2}}$$";
                specificCase.value="$$\chi="+X+"$$";
            else
                X=0.877/(esbeltezReduzida.value^2);
                specificCase.compareCase=2;
                specificCase.compare="$$\gamma_{0}>1.5$$";
                specificCase.equationString="$$\chi=\frac{0.877}{\gamma_{0}^{2}}$$";
                specificCase.value="$$\chi="+X+"$$";
                
            end
            
            fprintf(Standards.designReportFile,'%s',Standards.designReport.QsiSection(specificCase,esbeltezReduzida));
        end
        
        function Nrd=CompressSolve(Standards,X,Q,Ag,Fy,i)
           Nrd =(X*Q*Ag*Fy)/1.1;
           nrd.value=Nrd;
           nrd.equationValue="$$Nrd=\frac{"+X+"*"+Q+"*"+Ag+"*"+Fy+"}{"+1.1+"}$$"
           nrd.equationString="$$Nrd=\frac{X*Q*Ag*f_{y}}{\gamma_{a}}$$"
           fprintf(Standards.designReportFile,'%s',Standards.designReport.NrdSection(nrd));
           
            allStructuralGroups=getappdata(0,'allStructuralGroups');
            
             verifyInfo.type='Nrd';
               verifyInfo.rd=Nrd;
                allStructuralGroups(i).compressionVerify=verifyInfo;
                 %UPDATE STRUCTURALGROUPS AND MODEL
                                    setappdata(0,'allStructuralGroups',allStructuralGroups);
                                    model = getappdata(0,'model');
                                    model.allStructuralGroups=allStructuralGroups;
                                    setappdata(0,'model',model);
                                   
                          %UPDATE STRUCTURALGROUPS AND MODEL

        end
        
        
        
       
        %%CORTANTE
          function [lambda,lambdaP,lambdaR]=shearResistence(Standards,element,geometricType,designReportFile,i)
              
     
       
            switch geometricType
                case 1
                 
                    lambda.value=(element.section.h)/element.section.tw;
                    lambdaP.value=1.10*sqrt(5*element.material.elasticity/element.material.leakage);
                    lambdaR.value=1.37*sqrt(5*element.material.elasticity/element.material.leakage);
                    vpl.value=0.60*element.section.d*element.section.tw*element.material.leakage;
                    
                     lambda.equationString="$$\lambda=\frac{b}{t}$$";
                     lambdaP.equationString="$$\lambda_p=1.10*\sqrt{\frac{k_v*E}{f_y}}$$";
                     lambdaR.equationString="$$\lambda_r=1.37*\sqrt{\frac{k_v*E}{f_y}}$$";
                     lambda.equationValue="$$\lambda=\frac{"+element.section.h+"}{"+element.section.tw+"}$$";
                     lambdaP.equationValue="$$\lambda_p=1.10*\sqrt{\frac{"+5+"*"+element.material.elasticity+"}{"+element.material.leakage+"}}$$";
                     lambdaR.equationValue="$$\lambda_r=1.37*\sqrt{\frac{"+5+"*"+element.material.elasticity+"}{"+element.material.leakage+"}}$$";
                     
                     vpl.equationString="$$Vpl=0.60*Aw*F_y$$";
                     vpl.equationValue="$$Vpl=0.60*"+element.section.d*element.section.tw+"*"+element.material.leakage+"$$";
                    
                     
                     if( lambda.value<=  lambdaP.value)
                         specificCase.compareCase=1;
                         specificCase.compareString="$$\lambda\leq \lambda _{p}$$";
                         specificCase.compareValue="$$"+lambda.value+"\leq "+lambdaP.value+"$$";
                         vrd.value=vpl.value/1.1;
                         vrd.equationString="$$V_{rd}=\frac{V_{pl}}{\gamma _{a1}}$$";
                         vrd.equationValue="$$V_{rd}=\frac{"+vpl.value+"}{"+1.1+"}$$";
                     elseif (lambdaP.value<lambda.value && lambda.value <= lambdaR.value )
                         specificCase.compareCase=2;
                         specificCase.compareString="$$\lambda_{p}< \lambda \leq \lambda _{r}$$";
                         specificCase.compareValue="$$\"+lambdaP.value+"< \"+lambda.value+" \leq \"+lambdaR.value+"$$";
                         vrd.value=(lambdaP.value/lambda.value)*(vpl.value/1.1);
                         vrd.equationString="$$V_{rd}=\frac{\lambda_{p}}{\lambda }*\frac{V_{pl}}{\gamma_{a1}}$$";
                         vrd.equationValue="$$V_{rd}=\frac{"+lambdaP.value+"}{\"+lambda.value+" }*\frac{"+vpl.value+"}{"+1.1+"}$$";
                     else
                         specificCase.compareCase=3;
                         specificCase.compareString="$$\lambda >\lambda _{r}$$";
                         specificCase.compareValue="$$\"+lambda.value+" >"+lambdaR.value+"$$";
                         vrd.value=1.24*((lambdaP.value/lambda.value)^2)*(vpl.value/1.1);
                         vrd.equation="$$V_{rd}=1.24*(\frac{\lambda_{p}}{\lambda})^{2}*\frac{V_{pl}}{\gamma_{a1}}$$";
                         vrd.equationValue="$$V_{rd}=1.24*(\frac{\"+lambdaP.value+"}{\"+lambda.value+"})^{2}*\frac{"+vpl.value+"}{"+1.1+"}$$";
                     end
            
                        
%                       if(Vpl>=element.maxShearForce_XY*1000)
%                         data.approves=true;
%                         disp("aqui")
%                       else
%                         disp("aqui2")
%                         data.approves=false;
%                       end
%             


                        verifyInfo.type='Vrd';
                        verifyInfo.rd=vrd.value;
                         allStructuralGroups=getappdata(0,'allStructuralGroups');
                       
                         if (vrd.value>=(allStructuralGroups(i).maxShearForce_XY*1000))
                             vrd.approves=true;
                             
                         else
                             vrd.approves=false;
                             
                         end
                        allStructuralGroups(i).shearXYVerify=verifyInfo;
                        
                                    %UPDATE STRUCTURALGROUPS AND MODEL
                                    setappdata(0,'allStructuralGroups',allStructuralGroups);
                                    model = getappdata(0,'model');
                                    model.allStructuralGroups=allStructuralGroups;
                                    setappdata(0,'model',model);
                                   
                                    %UPDATE STRUCTURALGROUPS AND MODEL
                        
                     fprintf(designReportFile,'%s',Standards.designReport.VrdSection(lambda,lambdaP,lambdaR,vpl,vrd,specificCase));
                    
             
                    
                otherwise
            end
         
            
            
        end
        
        
        %%FUNÇÃO PARA VERIFICAÇAO TRAÇÃO
        function NrdSolve(Standards,element,designReportFile,i)
        allStructuralGroups=getappdata(0,'allStructuralGroups');
            switch element.section.Name
                case "WBeam"
                    NrdB=((element.section.Ag)*((element.material.leakage)))/1.1;%m2,n/m2
                    NrdB=NrdB*1e-3;%transform in Kn
                    
                    if NrdB>element.maxTraction
                        approves=true;
                        
                    else
                        approves=false;
                        
                    end
                    
                    if element.maxTraction<1e-10
                        element.maxTraction=0;
                    
                        
                    end
                    Nrd.rd=NrdB;
                    Nrd.type="Nrdb";
                  allStructuralGroups(i).tractionVerify=Nrd;
     %UPDATE STRUCTURALGROUPS AND MODEL
            setappdata(0,'allStructuralGroups',allStructuralGroups);
            model = getappdata(0,'model');
            model.allStructuralGroups=allStructuralGroups;
            setappdata(0,'model',model);
             fprintf(designReportFile,'%s', Standards.designReport.tractionSection(1.1,i,Nrd.rd,element.section.Ag,element.material.leakage,approves,element.maxTraction));
                 
   %UPDATE STRUCTURALGROUPS AND MODEL
                    
            end
        end
        
        
        function axialBendingCombination(~,element)
%            allStructuralGroups=getappdata(0,'allStructuralGroups');
%            
%          teste=[];
%        
%        
%           for z=1:size(element.bendingXYVerify,2)
%              tractionCombination=[];
%              for i=1:size(element.tractionCases,2)
%             
%                  for j=1:size(element.tractionCases(i).result,2)
%                       
%                     if element.tractionCases(i).result(j)/element.tractionVerify.rd>=0.2
%                         verification=(element.tractionCases(i).result(j)/element.tractionVerify.rd)+(8/9)*(abs(element.bendingCases_XY(i).result(j)*1000)/element.bendingXYVerify(z).rd)+0;
%                         tractionCombination(i,j)=verification;
%                     else
% 
%                       verification=((element.tractionCases(i).result(j))/(2*element.tractionVerify.rd))+(abs((element.bendingCases_XY(i).result(j)*1000))/(element.bendingXYVerify(z).rd))+0;
%                        tractionCombination(i,j)=verification;
%                     end
%                  end
%              end 
%              teste=[teste;tractionCombination];
%           end
%           
%            teste
          
        end
        
        
       
        
        
        %         %%FUNCÕES FLEXÃO SIMPLES
        %
        % EIXO DE MAIOR INERCIA DOIS EIXOS DE SIMETRIA PARA ALMA NAO
        % ESBELTA
        
        function flt = fltSolver(Standards,element,fid,i) %Calcula a Flambagem lateral com torção, no eixo forte em [N.mm] SEÇÃO G.2.1
            b1 = 0.7 * element.section.Wx*element.material.leakage/((element.material.elasticity * 1.72e-8)); %PAG135  element.section.J subs por 1.72x10-8CONFERIR FORMULA ADE J e B1
            
            %[M,~] = element.intBendingMoment_XY([(element.length/4),(element.length/2),((element.length*3)/4)]);
            
            position1=normalizeLength(element,4);
            position2=normalizeLength(element,2);
            position3=normalizeLength(element,(3/4));
        
            [M,~] = element(position1(1)).intBendingMoment_XY([(position1(2)),(position2(2)),(position3(2))]);
            disp('CB AQUI');
            cb = element.maxBendMoment_XY(1) *12.5 /(element.maxBendMoment_XY(1)*2.5 + M(1)*3 + M(2)*4 + M(3)*3)%PERGUNTAR PARA O PEDRO A EXTREMIDADE QUE COMEÇA A DIVISAO e com pegar exato pag 47 nbr
            %cb=maxmoment*12.5/maxmoment*2.5+momentu1/4 lenght*3+momentu2/4
            %lenght*3+momentu3/4 lenght*3  apartir da esquerda
            if cb > 3
                cb = 3;
            end
            disp("-----------------------------------LATERAL------------------------------------------")
            lambidaTorcao = element.length/element.section.ry %element.lenght PROVISORIO trocar por LB
            %lambidaTorcao=distancia do travemnto lateral/raio de giraçao y
            lambidaLimite1 = 1.76*sqrt(element.material.elasticity/element.material.leakage)
            %lambidaLimite1 = 1.76*sqrt(Gparam.Ea/Gparam.fy);
            lambidaLimite2 = 1.38*sqrt(element.section.Iy * 1.72e-8)/(element.section.ry * 1.72e-8 * b1) * sqrt( 1 + sqrt( 1 + (27*element.section.cw * b1^2 / element.section.Iy)))%element.section.J subs por 1.72x10-8CONFERIR FORMULA ADE J
            %lambidaLimite2 = 1.38*sqrt(perfil.Iy * perfil.J)/(perfil.ry * perfil.J * b1) * sqrt( 1 + sqrt( 1 + 27*perfil.cw * b1.^2 / perfil.Iy));
            if lambidaTorcao < lambidaLimite1
                %               flt = perfil.mpl;
                data.lambdaCompare=" $$\lambda < lambda p.$$ ";
                flt = element.section.mpl;
                data.mcrForm=" $$Mrd={Mpl \over \gamma a1}.$$";
                data.mcrFormSubs=" $$Mrd={" + element.section.mpl + " \over 1.1}.$$";
                
            elseif lambidaTorcao < lambidaLimite2
                data.lambdaCompare=" $$ \lambda p < \lambda < \lambda r .$$";
                %               flt = cb * (perfil.mpl - (perfil.mpl - 0.7*perfil.Wx*Gparam.fy)*(lambidaTorcao - lambidaLimite1)/(lambidaLimite2 - lambidaLimite1));
                flt = cb * (element.section.mpl - (element.section.mpl -  0.7*element.section.Wx*element.material.leakage)*(lambidaTorcao - lambidaLimite1)/(lambidaLimite2 - lambidaLimite1));
                data.McrForm="$$Mrd={{Cb\over \gamma a1} \star[{Mpl -(Mpl-Mr)\star {\lambda -\lambda p \over \lambda r - \lambda p} }]}.$$";
                data.McrFormSubs="$$Mrd={{"+ cb + "\over \1.1} \star[{" + element.section.mpl + " - " + element.section.mpl + "-" +(0.7*element.section.Wx*element.material.leakage)+ "\star \" + lambidaTorcao + " - " + lambidaLimite1 + " \over \" + lambidaLimite2 + " - " + lambidaLimite1 + "} }]}.$$";
                data.cb=cb;
            else
                %               flt = (cb * pi.^2 * Gparam.Ea * perfil.Iy) / Distparam.lb.^2 * sqrt( perfil.cw/perfil.Iy * (1 + 0.039 * perfil.J * Distparam.lb.^2/perfil.cw) );
                data.lambdaCompare=" $$\lambda > \lambda r .$$ ";
                flt = ((cb * pi^2 * element.material.elasticity * element.section.Iy) / element.length^2) * sqrt(element.section.cw/element.section.Iy * (1 + 0.039 * 1.721e-8 * element.length^2/element.section.cw) );%element.section.J subs por 1.72x10-8CONFERIR FORMULA ADE J
                
                data.mcrForm='$$Mcr={{Cb \star \pi^{2} \star E \star Iy\over Lb^{2}} \sqrt{{Cw\over Iy}(1+0.039{J\star Lb^{2} \over Cw})}}$$';
                data.mcrFormSubs="$$Mcr={{" + cb + " \star \pi^{2} \star " + element.material.elasticity + " \star " + element.section.Iy + "\star" + element.length + "^{2}} \sqrt{{" + element.section.cw + "\over Iy}(1+0.039{"+ 1.721e-8 +"\star" + element.length +"^{2} \over "+ element.section.cw + "})}}$$";
                
                
            end
          
            cb=cb
            mr=0.7*element.section.Wx*element.material.leakage
            disp("flt aqui");
           
            flt
            
            %data to report
            %             data.Lambda= "$$\lambda={ Lb\over ry 1.1}.$$";
            %             data.LambdaP= "$$Nrd={1.76 \star \sqrt{E \over Fy}}.$$";
            %             data.LambdaR=" $$Mcr={{1.38 \sqrt{IyJ} \over ry \star J \star \beta} \star \sqrt{1+\sqrt{1+27 \star Cw \star \beta^{2} \over Iy}}}$$";
            data.Lambda="$$\lambda={ "+ element.length +"\over " + element.section.ry + "}.$$";
            data.LambdaP="$$\lambda p={1.76 \star \sqrt{" +element.material.elasticity + " \over " +element.material.leakage + "}}.$$";
            data.LambdaR="$$\lambda r={{1.38 \sqrt{" +element.section.Iy +" \star"+ 1.721e-8 +"} \over " + element.section.ry +" \star "+ 1.721e-8 +" \star \beta} \star \sqrt{1+\sqrt{1+27 \star " + element.section.cw + " \star \beta^{2} \over " + element.section.Iy + "}}}$$";
            data.bl=b1;
            data.M=M;
            data.lambidaTorcao=lambidaTorcao;
            data.lambidaLimite1=lambidaLimite1;
            data.lambidaLimite2=lambidaLimite2;
            data.flt=flt;
            %             data.element=element;
            data.approves=true;
            
            
            %             fprintf(fid,'%s', Standards.designReport.flexSection(i,data));
        end
        
        function flm = flmSolver(Standards,element) %Calcula a Flambagem lateral de mesa, no eixo forte em [N.mm] SEÇÃO G.2.2
            disp("-----------------------------------MESA------------------------------------------")
            lambidaMesa = (element.section.bf/2)/element.section.tf
            lambidaLimite1 = 0.38*sqrt(element.material.elasticity/element.material.leakage)
            %           lambidaLimite1 = 0.38*sqrt(Gparam.Ea/Gparam.fy);
            lambidaLimite2 = 0.83*sqrt(element.material.elasticity/(0.7*element.material.leakage))%TABELA G.1 >>NOTA 6
            %           lambidaLimite2 = 0.83*sqrt(Gparam.Ea/(0.7*Gparam.fy));
            if lambidaMesa < lambidaLimite1
                flm = element.section.mpl;
            elseif lambidaMesa < lambidaLimite2
                flm = element.section.mpl - (element.section.mpl - 0.7*element.section.Wx*element.material.leakage)*(lambidaMesa - lambidaLimite1)/(lambidaLimite2 - lambidaLimite1);
            else
                flm = (0.69 * element.material.elasticity * element.section.Wx) / lambidaMesa.^2;
            end
            disp("flm="+flm);
        end
        
        function fla = flaSolver (Standards,element) %Calcula a Flambagem lateral de Alma, no eixo forte em [N.mm] SEÇÃO G.2.2
            disp("-----------------------------------ALMA------------------------------------------")
            lambidaAlma = element.section.h/element.section.tw
            %             lambidaLimite1 = 3.76*sqrt(Gparam.Ea/Gparam.fy);
            lambidaLimite1 = 3.76*sqrt(element.material.elasticity/element.material.leakage)
            %             lambidaLimite2 = 5.70*sqrt(Gparam.Ea /Gparam.fy);
            lambidaLimite2 = 5.70*sqrt(element.material.elasticity /element.material.leakage)
            if lambidaAlma < lambidaLimite1
                fla = element.section.mpl;
            elseif lambidaAlma < lambidaLimite2
                fla = element.section.mpl - (element.section.mpl - 0.7*element.section.Wx*element.material.leakage)*(lambidaAlma - lambidaLimite1)/(lambidaLimite2 - lambidaLimite1);
            else
                element.section.AlmaEsbelta = 1;
                fla = 1000000000000000000;
            end
            disp("fla="+fla);
        end
        
        
        %EIXO DE MAIOR INERCIA DOIS EIXOS DE SIMETRIA PARA CASO ALMA
        % ESBELTA
        function flm = flmEsbeltSolver (element) %Calcula a Flambagem lateral da mesa, no eixo forte caso a alma seja esbelta em [N.mm]
            %                At = element.h*perfil.tw/(perfil.tf*perfil.bf);
            At = element.section.h*element.section.tw/(element.section.tf*element.section.bf);
            if At > 10
                set (handles.text1, 'String', 'Perfill W com alma Extremamente Esbelta');
                return
            else
                %                   kpg = 1 - At/(1200+300*At)*(perfil.h/perfil.tw - 5.70*sqrt(Gparam.Ea/Gparam.fy));
                kpg = 1 - At/(1200+300*At)*(element.section.h/element.section.tw - 5.70*sqrt(element.material.elasticity/element.material.leakage));
                %                   lambidaMesa = (perfil.bf/2)/perfil.tf;
                lambidaMesa = (element.section.bf/2)/element.section.tf;
                %                   lambidaLimite1 = 0.38*sqrt(Gparam.Ea/Gparam.fy);
                lambidaLimite1 = 0.38*sqrt(element.material.elasticity/element.material.leakage);
                %                   kc = 4/sqrt(perfil.h/perfil.tw);
                kc = 4/sqrt(element.section.h/element.section.tw);
                if kc < 0.35
                    kc = 0.35;
                elseif kc > 0.76
                    kc = 0.76;
                end
                %                   lambidaLimite2 = 0.95*sqrt(Gparam.Ea * kc /(0.7*Gparam.fy));
                lambidaLimite2 = 0.95*sqrt(element.material.elasticity * kc /(0.7*element.material.leakage));
                if lambidaMesa < lambidaLimite1
                    %                        flm = kpg * perfil.Wx * Gparam*fy;
                    flm = kpg * element.section.Wx * element.material.leakage;
                elseif lambidaMesa < lambidaLimite2
                    %                         flm = (kpg * perfil.Wx * Gparam.fy) * (1 - 0.3 *(lambidaMesa - lambidaLimite1)/(lambidaLimite2 - lambidaLimite1));
                    flm = (kpg * element.section.Wx * element.material.leakage) * (1 - 0.3 *(lambidaMesa - lambidaLimite1)/(lambidaLimite2 - lambidaLimite1));
                else
                    %                         flm = (kpg * 0.90 * Gparam.Ea * kc * perfil.Wx) / lambidaMesa.^2;
                    flm = (kpg * 0.90 * element.material.elasticity * kc * element.section.Wx) / lambidaMesa.^2;
                end
            end
        end
        
        function flt = fltEsbeltSolver (Standards,element) %Calcula a Flambagem lateral com torção, no eixo forte caso a alma seja esbelta em [N.mm]
            
            [M,~] = element.intBendingMoment_XY([(element.length/4),(element.length/2),(element.length/3)]);
            cb = element.maxBendMoment_XY *12.5 /(element.maxBendMoment_XY*2.5 + M(1)*3 + M(2)*4 + M(3)*3);
            
            if cb > 3
                cb = 3;
            end
            %            At = perfil.h*perfil.tw/(perfil.tf*perfil.bf);
            At = element.section.h*element.section.tw/(element.section.tf*element.section.bf);
            if At > 10
                set (handles.text1, 'String', 'Perfill W com alma Extremamente Esbelta');
                return
            else
                %                   kpg = 1 - At/(1200+300*At)*(perfil.h/perfil.tw - 5.70*sqrt(Gparam.Ea/Gparam.fy));
                kpg = 1 - At/(1200+300*At)*(element.section.h/element.section.tw - 5.70*sqrt(element.material.elasticity/element.material.leakage));
                %                   Iy2 = (perfil.tf*perfil.bf.^3)/12 + (perfil.h/3 * perfil.tw.^3)/12;
                Iy2 = (element.section.tf*element.section.bf.^3)/12 + (element.section.h/3 * element.section.tw.^3)/12;
                %                   Ag2 = (perfil.bf*perfil.tf + (perfil.h/3*perfil.tw));
                Ag2 = (element.section.bf*element.section.tf + (element.section.h/3*element.section.tw));
                ry2 = sqrt(Iy2/Ag2);
                %                   lambidaTorcao = Distparam.lb/ry2;
                lambidaTorcao = element.length/ry2
                %                   lambidaLimite1 = 1.10*sqrt(Gparam.Ea/Gparam.fy);
                lambidaLimite1 = 1.10*sqrt(element.material.elasticity/element.material.leakage)
                %                   lambidaLimite2 = pi * sqrt(Gparam.Ea/0.7*Gparam.fy);
                lambidaLimite2 = pi * sqrt(element.material.elasticity/0.7*element.material.leakage)
                if lambidaTorcao < lambidaLimite1
                    %                       flt = kpg * perfil.Wx * Gparam*fy;
                    flt = kpg * element.section.Wx * element.material.leakage;
                elseif lambidaTorcao < lambidaLimite2
                    %                       flt = (cb * kpg * perfil.Wx * Gparam.fy) * (1 - 0.3 *(lambidaTorcao - lambidaLimite1)/(lambidaLimite2 - lambidaLimite1));
                    flt = (cb * kpg * element.section.Wx * element.material.leakage) * (1 - 0.3 *(lambidaTorcao - lambidaLimite1)/(lambidaLimite2 - lambidaLimite1));
                else
                    %                         flt = (cb * kpg * pi.^2 * Gparam.Ea * perfil.Wx)/(lambidaTorcao.^2);
                    flt = (cb * kpg * pi.^2 * element.material.elasticity * element.section.Wx)/(lambidaTorcao.^2);
                end
                if flt > kpg * element.section.Wx * element.material.leakage      %flt > kpg * perfil.Wx * Gparam.fy
                    flt = kpg * element.section.Wx * element.material.leakage; %flt = kpg * perfil.Wx * Gparam.fy;
                end
            end
            
        end
        
        
        
        % EIXO DE MENOR INERCIA DOIS EIXOS DE SIMETRIA FLA NAO SE APLICA
        % PARA O EIXO DE MENOR INERCIA NOTA 3
        
        function flmy = flmySolver(element) %Calcula a Flambagem lateral de mesa, no eixo fraco em [N.mm]  SEÇÃO G.2.2
            %             lambidaMesa = (perfil.bf/2)/perfil.tf;
            lambidaMesa = (element.section.bf/2)/element.section.tf;
            
            %            lambidaLimite1 = 0.38*sqrt(Gparam.Ea/Gparam.fy);
            lambidaLimite1 = 0.38*sqrt(element.material.elasticity/element.material.leakage);
            lambidaLimite2 = 0.83*sqrt(element.material.elasticity/(0.7*element.material.leakage));
            if lambidaMesa < lambidaLimite1
                %               flmy = perfil.mply;
                flmy = element.section.mply;
            elseif lambidaMesa < lambidaLimite2
                %              flmy = perfil.mply - (perfil.mply - perfil.Wy*Gparam.fy)*(lambidaMesa - lambidaLimite1)/(lambidaLimite2 - lambidaLimite1);
                flmy = element.section.mply - (element.section.mply - element.section.Wy*element.material.leakage)*(lambidaMesa - lambidaLimite1)/(lambidaLimite2 - lambidaLimite1);
            else
                %              flmy = (0.69 * Gparam.Ea * perfil.Wy) / lambidaMesa.^2; %NOTA 6
                flmy = (0.69 * element.material.elasticity * element.section.Wy) / lambidaMesa.^2;
            end
        end
        
        
    end
    
end

