classdef loadAndUpdateInfos<handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
      handlesElementsLoads
      listagemDePerfil
      listagemDeCasosDeCarga
    
     
    end
    
       
   
    
    methods
        function this = loadAndUpdateInfos()
               this.handlesElementsLoads=guidata(findobj('Tag','GUI_ElementLoads'));
        end
     %REFERENTE A PERFIL%   
%         function loadStringInPerfilList(this,v)%v numero da linha do grupo de familias.
%                         %x=this.xsacomin{1}{1}{1,1}
%                        
%                         this.listagemDePerfil=[];
%                         for i=1:1:length(this.xsacomin{v})%for que caminha no numero de familias
%                             linhas=size(this.xsacomin{v}{i});%conta quantos perfis estao na familia corrente do for
%                             for j=1:1:linhas(1,1)%for que anda nas familias dento da familia corrente (i)
%                              
%                              this.listagemDePerfil=string([this.listagemDePerfil; this.xsacomin{v}{i}{j,1}]);%concatena o nome dos perfis de todas familias.
%                              
%                             end
%                         end
%                                              
%                    data1=cellstr(this.listagemDePerfil);%tranforma o vetor de string gerado em cell
%                    set(this.handles.listOfPerfil,'String',data1);%seta o valor fa list box
%                         
%                     switch v
%                             case 1
%                              set(this.handles.imgShape,'CData',imread('iBeam.jpg'));
%                     end
%         end
%          function loadPerfilProperties(this,v)
%                      data1=cell(4,2); 
%                         [j,i]=this.teste(v);
%                         data1(1,:)={'BF [mm]',this.xsacomin{v}{i}{j,3}};
%                         data1(2,:)={'TF [mm]',this.xsacomin{v}{i}{j,5}};
%                         data1(3,:)={'TW [mm]',this.xsacomin{v}{i}{j,4}};
%                         data1(4,:)={'Comprimento [m]','0'};
%                         this.handles.UITable.Data=data1;
%                         set(this.handles.uitable1,'Data',data1);
%                         
%                         
%          end
%          
%          function loadCarreamento(this)
%                      data1=cell(2,2); 
%                        
%                         data1(1,:)={'QX-Kn/m',0};
%                         data1(2,:)={'QY-Kn/m',0};
%                         
%                         this.handles.UITable.Data=data1;
%                         set(this.handles.uitable1,'Data',data1);
%                         
%                         
%          end
%          
%          function [index,i]=teste(this,v)
%             
%               index = get(this.handles.listOfPerfil,"Value");
%               a=0;
%               b=0;
%            
%               for i=1:1:size(this.xsacomin{v},2)
%                   
%                 a=a+size(this.xsacomin{v}{i},1)  ;
%                 if(a>=index)
%                     index=index-b;
%                     return
%                 else
%                     b=a;
%                 end
%               end
%         end
       %REFERENTE A PERFIL%   
       function setLoadCasesList(this,name,qx,qy,qz,qx1,qx2,qy1,qy2,qz1,qz2,dtx,dty,dtz,type,combination,class,psi)
           this.handlesElementsLoads=guidata(findobj('Tag','GUI_ElementLoads'));
           NewRow = { name,qx,qy,qz,qx1,qx2,qy1,qy2,qz1,qz2,dtx,dty,dtz,type,combination,class,psi};
           uimulticollist ( this.handlesElementsLoads.loadDefinedCasesList, 'addRow', NewRow );
           %set(this.handlesElementsLoads.loadNameCaseList,'String',NewRow(1))
           setappdata(0,'listdata',get(this.handlesElementsLoads.loadDefinedCasesList,'ApplicationData'));
       end
       
       function deleteRow(this,value)
           disp(value);
           uimulticollist ( this.handlesElementsLoads.loadDefinedCasesList, 'delRow', value);
           setappdata(0,'listdata',get(this.handlesElementsLoads.loadDefinedCasesList,'ApplicationData'));
           
       end
       
       function editRow(this,value,str)
           this.handlesElementsLoads=guidata(findobj('Tag','GUI_ElementLoads'));
           listdata=getappdata(0,'listdata');
           set(this.handlesElementsLoads.edit_Qx,'String',listdata.string(value,2));
           set(this.handlesElementsLoads.edit_Qy,'String',listdata.string(value,3));
           set(this.handlesElementsLoads.edit_Qz,'String',listdata.string(value,4));
           set(this.handlesElementsLoads.edit_Qx1,'String',listdata.string(value,5));
           set(this.handlesElementsLoads.edit_Qx2,'String',listdata.string(value,6));
           set(this.handlesElementsLoads.edit_Qy1,'String',listdata.string(value,7));
           set(this.handlesElementsLoads.edit_Qy2,'String',listdata.string(value,8));
           set(this.handlesElementsLoads.edit_Qz1,'String',listdata.string(value,9));
           set(this.handlesElementsLoads.edit_Qz2,'String',listdata.string(value,10));
           set(this.handlesElementsLoads.edit_dtx,'String',listdata.string(value,11));
           set(this.handlesElementsLoads.edit_dty,'String',listdata.string(value,12));
           set(this.handlesElementsLoads.edit_dtz,'String',listdata.string(value,13));
            length(str)
          listdata.string(value,15)
%            set(this.handlesElementsLoads.loadCombination,'String',str{set(this.handlesElementsLoads.loadCombination,'Value',this.check(str,listdata.string(value,15)))});%onde esta o 2 index da string 
         
       end
       
       
       function changePopUpOptions(this,value)
           variableLoadClassString={'Efeito Temperatura' 'Ação Vento' 'Acões Truncadas' 'Outros*'};
           permanentLoadClassString={'Peso própio de estr.metálica Peso própio de estr.pré-moldada' 'Peso própio de estr.moldadas no local/elem.construtivos industrializados e empuxos permanentes'...
               'Peso própio de elem.contrutivos industrializados com adições in loco' 'Peso própio de elem.construtivos  em geral e equipamentos' 'Indiretas'};
           switch value
               case 1
                   set(this.handlesElementsLoads.psi,'Enable','off');
                   set(this.handlesElementsLoads.psi,'String','Não Possui');
                   set(this.handlesElementsLoads.loadClass,'String',permanentLoadClassString);
               case 2
                   set(this.handlesElementsLoads.loadClass,'String',variableLoadClassString);
                   set(this.handlesElementsLoads.psi,'Enable','on');
                   set(this.handlesElementsLoads.psi,'String',{'teste1' 'teste2'});
           end
       end
       
      
        function i=check(str,stringSaved)
           disp("foi");
             for i=1:1:length(str)
                 if strcmp(stringSaved,str(i))
                     return 
                 end
             end
        end
       
    end
    
    
    
end

