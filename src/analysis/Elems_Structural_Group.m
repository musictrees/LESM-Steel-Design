classdef Elems_Structural_Group<handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name;
        elems=[];
        ids=[]
        length=0;
        Kx;
        Ky;
        Lx;
        Ly;
        Cb;
        Lb;
        Kz;
        
        maxShearForce_XY= []
        maxShearForce_XY_CombinationId=[];
        maxBendMoment_XY_CombinationId=[];
        maxBendMoment_XY= [];
        
        
        maxCompression=[];
        maxCompressionId=[];
        compressionVerify=[];
        
        
        maxTraction=[];
        maxTractionId=[];
        tractionVerify=[];
        
 
        
        
        maxShearForce_XZ= []
        
       
        
        
        maxBendMoment_XZ= []
        
       
        section;
        material;
        nodes;
        %% tracao e compressao
        axialForce=[];
        axialForceCases=[];
        %%
        %%sepearando so tracao por caso
        tractionCases=[];
        %%seprando comrpessao por caso
        compressionCases=[];
        
        %%
        bending_XY;% momentos salvos para cada tipo de carga
        bendingCases_XY;% momentos salvos para cada combinação
        bendingXYVerify; %guardo os fatores de dimensionamento
        
        shear_XY;% momentos salvos para cada tipo de carga
        shearCases_XY;% momentos salvos para cada combinação
        shearXYVerify;%guardo os fatores de dimensionamento
        
    end
    
    methods
        
        function this=Elems_Structural_Group(elems,name,ids,Kx,Ky,Lx,Ly,Cb,Lb,Kz)
            this.Kx=Kx;
            this.Ky=Ky;
            this.Lx=str2double(Lx);
            this.Ly=str2double(Ly);
            this.Cb=Cb;
            this.Lb=str2double(Lb);
            this.Kz=str2double(Kz);
            
            this.elems=elems;
            this.ids=ids;
            this.name=name;
            this.section=elems(1).section;%por enquanto assume a msm matersecçao pára todos
            this.material=elems(1).material;%por enquanto assume o msm material para todos
            this.transformElems();
            
        end
        %%traction
        function separeTractionCases(this)
           
            for i=1:size(this.axialForceCases,2)
                x=this.axialForceCases(i).result;
                
                x(1,x(1,:)<0)=0;
                CombinationId=this.axialForceCases(i).CombinationId;
                CombinationSubType=this.axialForceCases(i).CombinationSubType;
                
                this.tractionCases(i).result=x(1,:);
                this.tractionCases(i).CombinationId=CombinationId;
                this.tractionCases(i).CombinationSubType=CombinationSubType;
            end
        end
        
        function getMaxTraction(this)
              values=[];
            for i=1:1:size(this.tractionCases,2)
                values=[values;max(abs(this.tractionCases(i).result))];
            end
            this.maxTraction=max(values);
            combIdIndex=values==this.maxTraction;
            this.maxTractionId=this.tractionCases(combIdIndex).CombinationId;   
        end
        %%
          %%compress
        function separeCompressCases(this)
            
               for i=1:size(this.axialForceCases,2)
               
                 
                y=this.axialForceCases(i).result;
                CombinationId=this.axialForceCases(i).CombinationId
                CombinationSubType=this.axialForceCases(i).CombinationSubType;
                y(1,y(1,:)>0)=0; %pega valores menores q 0
               
                this.compressionCases(i).result=y(1,:);
                this.compressionCases(i).CombinationId=CombinationId;
                this.compressionCases(i).CombinationSubType=CombinationSubType;
   
              end
        end
        
        function getMaxCompress(this)
              values=[];
            for i=1:1:size(this.compressionCases,2)
                values=[values;max(abs(this.compressionCases(i).result))];
            end
            this.maxCompression=max(values);
            combIdIndex=values==this.maxCompression;
            this.maxCompressionId=this.compressionCases(combIdIndex).CombinationId;  
            
        end
        %%
        function getMaxBending(this)
           
            values=[];
            for i=1:1:size(this.bendingCases_XY,2)
                values=[values;max(abs(this.bendingCases_XY(i).result))];
            end
            this.maxBendMoment_XY=max(values);
            combIdIndex=values==this.maxBendMoment_XY;
            this.maxBendMoment_XY_CombinationId=this.bendingCases_XY(combIdIndex).CombinationId;            
        end
        
            function getMaxShear(this)
           
            values=[];
            for i=1:1:size(this.shearCases_XY,2)
                values=[values;max(abs(this.shearCases_XY(i).result))];
            end
            this.maxShearForce_XY=max(values);
            combIdIndex=values==this.maxShearForce_XY;
            this.maxShearForce_XY_CombinationId=this.shearCases_XY(combIdIndex).CombinationId;            
        end
        function transformElems(this)
            
            for i=1:length(this.elems)
                this.length=this.length+this.elems(i).length; %soma o tamanho de todos eloementos que compoem o grupo
                %%this.maxBendMoment_XY=cat(2,this.maxBendMoment_XY,this.getbending_moment_Z_list(this.elems(i))); %concatena todos os momentos do plano XY dos elementos chama a função pois as vezes os arrays sao vazios precisando de uma verificação
                this.maxShearForce_XY=cat(2,this.maxShearForce_XY,this.getshear_force_Y_list(this.elems(i)));%concatena todos os cortantes do plano XY dos elementos chama a função pois as vezes os arrays sao vazios precisando de uma verificação
               % this.allAxialForce = [this.allAxialForce  this.elems(i).axial_force(1)];
                
            end
            
            
            this.nodes=this.getNodes(this.elems);
            
            disp("aquiM"+this.maxBendMoment_XY);
            disp("aquiV"+this.maxShearForce_XY);
            this.maxBendMoment_XY=max(abs(this.maxBendMoment_XY)); %pega o maior valor de moemento na lista de momentos em modulo , arruam rpara monosimetricos pis o max negativo e positivo influenciam
            this.maxShearForce_XY=max(abs(this.maxShearForce_XY)); %pega o maior valor de cortante na lista de contante em modulo , arruam rpara monosimetricos pis o max negativo e positivo influenciam
            
            %this.maxCompress=this.allAxialForce(this.allAxialForce>=0);
            %this.maxTraction=this.allAxialForce(this.allAxialForce<=0);
            
            
        end
        
        function array=getbending_moment_Z_list(~,elem)
            array=0;
            
            if(~isempty(elem.maxBendMoment_XY))
                array=cat(2,array,elem.maxBendMoment_XY(1));
                disp("entrou");
            else
                
            end
            
            if(~isempty(elem.bending_moment_Z))
                array=cat(2,array,elem.bending_moment_Z(1),elem.bending_moment_Z(2));%case do you wanna a correct direction multiply the first node per -1
                disp("entrou");
            else
                
            end
        end
        
        function array=getshear_force_Y_list(~,elem)
            array=0;
            
            if(~isempty(elem.maxShearForce_XY))
                array=cat(2,array,elem.maxShearForce_XY(1));
                disp("entrou");
            else
                
            end
            
            if(~isempty(elem.shear_force_Y))
                array=cat(2,array,elem.shear_force_Y(1),elem.shear_force_Y(2));
                disp("entrou");
            else
                
            end
        end
        %          function array=getshear_axial_force_list(~,elem)
        %             array=0;
        %
        %             if(~isempty(elem.maxAxialForce))
        %                 array=cat(2,array,elem.maxAxialForce(1));
        %                 disp("entrou");
        %             else
        %
        %             end
        %
        %             if(~isempty(elem.axial_force))
        %                 array=cat(2,array,elem.axial_force(1),elem.axial_force(2));
        %                 disp("entrou");
        %             else
        %
        %             end
        %         end
        function nodes=getNodes(~,elems)
            nodes=[elems(1).nodes(1),elems(length(elems)).nodes(2)]; %esa funçao pega o primeiro e o ultimo no do elemento total composto pelo grupo de elementos
        end
        
    end
end

