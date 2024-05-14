classdef DesignReport<handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods 
%         function fid=criarArquivo(~)
%             fid = fopen('teste_site.html','w');
%            
%           
%         end
        
        function html=criarHtml(~)
         
            html=["<!DOCTYPE html>"
            "<html lang='pt-br'>"

            "<head>"

              "<meta charset=”UTF-8”>"
              "<meta name='viewport' content='width=device-width, initial-scale=1, shrink-to-fit=no'>"
              "<meta name='description' content=''>"
              "<meta name='author' content=''>"

              "<title>Memorial de Calculo</title>"

              "<!-- Bootstrap core CSS -->"
              "<link href='vendor/bootstrap/css/bootstrap.min.css' rel='stylesheet'>"

              "<!-- Custom styles for this template -->"
              "<link href='css/simple-sidebar.css' rel='stylesheet'>"
              "<script src='vendor/jquery/jquery.min.js'></script>"
              "<script src='vendor/bootstrap/js/bootstrap.bundle.min.js'></script>"
              "<script src='https://polyfill.io/v3/polyfill.min.js?features=es6'></script>"
              "<script id='MathJax-script' async src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>"
              "<script src='js/myOwnScript.js'></script>"                
            "</head>"

            "<body>"
            
            "<div class='d-flex' id='wrapper'>"
            ];

            
        end
            
          function sideMenu = createSideMenu(this,n)
            sideMenu=["<div class='bg-light border-right' id='sidebar-wrapper'>"
                "<div class='sidebar-heading text-center'>Elementos </div>"
                "<div class='list-group list-group-flush'>"
                "<a href='#' class='list-group-item list-group-item-action bg-light' onclick='showLinks(this)'><i class='glyphicon glyphicon-chevron-right'></i>Tracao</a>"
                "<div id='elementosTracaoLinks' class='list-group list-group-flush'>" 
                 "<a href='#' class='list-group-item list-group-item-action bg-light'>elemento-1</a>"
                "<a href='#' class='list-group-item list-group-item-action bg-light'>elemento-2</a>"
                "<a href='#' class='list-group-item list-group-item-action bg-light'>elemento-3</a>"
                "</div>"
               
               
                
                "<a href='#' class='list-group-item list-group-item-action bg-light' onclick='showLinks(this)'><i class='glyphicon glyphicon-chevron-right'></i>Compressao</a>"
                "<div id='elementosCompressaoLinks' class='list-group list-group-flush'>"
                "<a href='#' class='list-group-item list-group-item-action bg-light'>elemento-1</a>"
                "<a href='#' class='list-group-item list-group-item-action bg-light'>elemento-2</a>"
                "<a href='#' class='list-group-item list-group-item-action bg-light'>elemento-3</a>"
                "</div>"
                
                
                "<a href='#' class='list-group-item list-group-item-action bg-light' onclick='showLinks(this)'><i class='glyphicon glyphicon-chevron-right'></i>Flexao</a>"
                "<div id='elementosFlexaoLinks' class='list-group list-group-flush'>"+ this.aLink(n) +  "</div>"         
              
                
                "<a href='#' class='list-group-item list-group-item-action bg-light' onclick='showLinks(this)'><i class='glyphicon glyphicon-chevron-right'></i>Cortante</a>"
                "<div id='elementosCortanteLinks' class='list-group list-group-flush'>"
                "<a href='#' class='list-group-item list-group-item-action bg-light'>elemento-1</a>"
                "<a href='#' class='list-group-item list-group-item-action bg-light'>elemento-2</a>"
                "<a href='#' class='list-group-item list-group-item-action bg-light'>elemento-3</a>"
                "</div>"
                "</div>"
                "</div>"];
            
        end
        

        
        function itens=aLink(~,n)
           a="";
           for i=1:1:n
            
             
                a=a+"<a href='#Flex-"+i+"' class='list-group-item list-group-item-action bg-light'>elemento-" + num2str(i) + "</a>";
              
           end
      itens=a;
    
        end
        
        function tractionDiv=tractionSection(~,n,i,NrdB,Ag,Fy,approves,Nsd)
           
            
            if approves==true
                approves= "<div class='alert alert-success text-center' role='alert'><h4 class='alert-heading'>Seção Aprovada</h4></div>";
                            
            else
                approves= "<div class='alert alert-danger text-center' role='alert'><h4 class='alert-heading'>Seção Reprovada</h4></div>";
            
            end
            
           tractionDiv=[
            "<div id='Traction-" + i + "' class='card border-dark col-md-12 bg-light ' style='padding: 0px;'>"
            "<div class='card-header font-weight-bold'>Tração</div>"
            "<div class='card-body'>"
            "<h5 class='card-title'>Elemento-" + i + "</h5>"
            "<p class='card-text'></p>Escoamento da seção Bruta:"
            "$$Nrd={A \star Fy\over \gamma ay}.$$"
            "</p>"
            "<p>"
            "$$Nrd={" + Ag*1e4 + " \star " + (Fy/1e7) + "\over \ 1.1}.$$"
            "</p>"
            "<p><em>"
            "A=Area bruta</br>"
            "Fy=Tensão de escoamento do aço</br>"
            "?au=Coef.Redução" + n + " "
            "</em></p>"
            
            "<p>"
            "$$Nrd>Nsd.$$"
            "</p>"
            
            "<p>"            + NrdB + ">" + Nsd +            "</p>"
            

            "<p class='card-text'></p>Escoamento da seção Líquida:"
            "$$Nrd={Ae \star Fy\over \gamma au}.$$</p>"  + approves + ...
            "</div></div>"];
           
         
            
        end
        %----------------------------------------------------------flexao simples
        function flexDiv=flexSection(~,lambda,lambdaP,lambdaR,Mr,Mcr,i,data)     
            
             
            
           flexDiv=[
            
             "<p class='card-text font-weight-bold'>"+data.title+"</p>" 
             "<h6>Tabela G.1-NBR8800-pág:144</h6>"
            
             "<div class='row border border-secondary'>" 
             "<div class='col-md-4'>" 
             "<p>" + data.Lambda + "</p>" 
             "</div>" 
             "<div class='col-md-4'>" 
             "<p>" 
             "<p>" + data.LambdaP + "</p>" 
             "</p>" 
             "</div>" 
             "<div class='col-md-4'>"
             "<p>"
             "<p>" + data.LambdaR + "</p>" 
             "</p>" 
             "</div>" 
             "</div>" 
             "<div class='row'>"
             "<div class='col-md-4'>" 
             "<p>" 
             "$$\lambda={" + lambda + "}$$" 
             "</p>" 
             "</div>" 
             "<div class='col-md-4'>" 
             "<p>" 
             "$$\lambda p={" + lambdaP + "}$$"  
             "</p>" 
             "</div>" 
             "<div class='col-md-4'>" 
             "<p>" 
             "$$\lambda r={" + lambdaR + "}$$" 
             "</p>" 
             "</div>" 
             "</div>"  
             ];    
        end
        function flexDiv=flexMrd(~,data)
            
            
                        
            if data.approves==true
                approves= [
                  
                    "<div class='alert alert-success text-center' role='alert'><h4 class='alert-heading'>Seção Aprovada</h4></div>"
                    ];
                            
            else
                approves= [
                  
                    "<div class='alert alert-danger text-center' role='alert'><h4 class='alert-heading'>Seção Reprovada</h4></div>"
                    ];
            
            end
%        
            flexDiv=[
             "<h6>Secão-G.2.2-pág:130</h6>"
             "<div class='row border border-secondary'>"
                 "<div class='col-md-12'>" 
                     "<p>" + data.MrdString +  "</p>" 
                 "</div>" 
             "</div>"
             "<div class='row'>"
                "<div class='col-md-12 text-center'>" 
                    "<p><span class='font-weight-bold'>Mrd</span>="+ data.Mrd + "</p>" 
                "</div>"
                   "<div class='col-md-12 text-center'>" 
                approves
                    "</div>"
              "</div>"
              
             ];
        end
        
        function htmlfechar=fecharHtml(~)
           htmlfechar=["</div></body></html>"];
        end
        
        
           %%CORTANTE
       
           
            
           function VrdDiv=VrdSection(~,lambda,lambdaP,lambdaR,vpl,vrd,specificCase)
               
               switch specificCase.compareCase
                   case 1
                       compare=[specificCase.compareString
                           specificCase.compareValue];
                         solve=[
                           vrd.equationString
                            vpl.equationString
                            vpl.equationValue
                            vrd.equationValue
                            vrd.value
                           ];
                   case 2
                      compare=[specificCase.compareString 
                          specificCase.compareValue];
                        solve=[
                           vrd.equationString
                            vpl.equationString
                            vpl.equationValue
                            vrd.equationValue
                            vrd.value
                           ];
                       
                   case 3
                       
                      compare=[specificCase.compareString
                          specificCase.compareValue];
                        solve=[
                           vrd.equationString
                            vpl.equationString
                            vpl.equationValue
                            vrd.equationValue
                            vrd.value
                           ];
                       
                   otherwise
               end
               
               
                if vrd.approves==true
                approves= "<div class='alert alert-success text-center' role='alert'><h4 class='alert-heading'>Seção Aprovada</h4></div>";
                            
            else
                approves= "<div class='alert alert-danger text-center' role='alert'><h4 class='alert-heading'>Seção Reprovada</h4></div>";
            
            end
%        
            aprovesDiv=[
           
           
                "<div class='col-md-12 text-center'>" 
                  approves
                  "</div>"
                 
              
             ];
            VrdDiv=[
                
            
            "<div class='col-md-12'>"
            "<p class='card-text font-weight-bold'>Fator Vrd-Y</p>"
            "<h6>Seção-5.4.3.1-NBR8800-pág:50</h6>"
            "</div>"
            "<div class='row border border-secondary'>"
            
            
            
            "<div class='col-md-4 text-center'>"
            "<p class='font-weight-bold'>Esbeltez</p>"
            "<p>"+ lambda.equationString+"</p>"
            
            "</div>"
            "<div class='col-md-4 text-center'>"
            "<p class='font-weight-bold'>Esbeltez Lim</p>"
            "<p>"+ lambdaP.equationString+"</p>"
            
            "</div>"
            "<div class='col-md-4 text-center'>"
            "<p class='font-weight-bold'>Esbeltez Lim</p>"
            "<p>"+ lambdaR.equationString+"</p>"
            
            "</div>"
            
            
            
            "<div class='col-md-4 text-center'>"
            
            "<p>"+ lambda.equationValue+"</p>"
            
            "</div>"
            "<div class='col-md-4 text-center'>"
            
            "<p>"+ lambdaP.equationValue+"</p>"
            
            "</div>"
            "<div class='col-md-4 text-center'>"
            
            "<p>"+ lambdaR.value+"</p>"
            
            "</div>"
             "<div class='col-md-4 text-center'>"
            
            "<p>"+ lambda.value+"</p>"
            
            "</div>"
            "<div class='col-md-4 text-center'>"
            
            "<p>"+ lambdaP.value+"</p>"
            
            "</div>"
            "<div class='col-md-4 text-center'>"
            
            "<p>"+ lambdaR.value+"</p>"
            
            "</div>"
            "<div class='col-md-12 text-center'>"
            "<p>Como:</p>"
            compare
            
            "<p>Logo:</p>"
            solve
            "</div>"
            aprovesDiv
            "</div>"
            
            
            ];
           end
        
        
         
            
       
        
        
        function openSectionHeader=openHeaderSection(~,i,title)
            openSectionHeader=[
            "<div id='Compres-" + i + "' class='card border-titledark col-md-12 bg-light mt-5' style='padding: 0px;'>" 
             "<div class='card-header font-weight-bold'>"+title+"</div>" 
              
              "<div class='row align-items-center  p-3'>" 
               "<div class='col-md-12'>" 
               "<h5 class='card-title'>Elemento-" + i + "</h5>" 
               "</div>"
            ];
        end
        
        function openCardDiv=openCardSection(~,i,title,maxStress,masStressId)
            
            openCardDiv=[ "<div id='Compres-" + i + "' class='card border-titledark col-md-12 bg-light mt-5' style='padding: 0px;'>" 
             "<div class='card-header font-weight-bold'>"+title+"</div>" 
         
              "<div class='row align-items-center  p-3'>" 
               "<div class='col-md-12'>" 
               "<h5 class='card-title'>Elemento-" + i + "</h5>"
               "<h6 class='card-title'>Caso de carga crítico:" + masStressId + "  logo valor critico:" + maxStress + "</h6>"
               
               "</div>"
               "<div class='card-body '>" ];
        end
        function openSectionDiv=openSection(~)
            openSectionDiv=[
              "<div class='card-body'>"
              ];
        end
        %%compress
        function kFactorDiv=kFactorSection(~,data)
         
%           switch  data.boundary
%               case "free/fixed"
%                   imgAdress="free_fixed.png";
%                   
%               otherwise
%                   imgAdress="none.png";
%           end
          
            
           kFactorDiv=[
           
             "<p class='card-text font-weight-bold'>Condições de Contorno</p>" 
             "<h6>Seção-5.4.3.1-NBR8800-pág:50</h6>"
            
             "<div class='row border border-secondary'>" 
             "<div class='col-md-6 text-center'>" 
             "<p>Direção-X</p>" 
            
             "</div>" 
             "<div class='col-md-6 text-center'>" 
             "<p>Direção-Y</p>" 
           
             "</div>" 
              
            
            
             "<div class='col-md-6 text-center'>" 
             "<p>" 
             "$$K_{x}={" + data.KxValue + "}$$" 
             "</p>" 
             "</div>" 
             "<div class='col-md-6 text-center'>" 
             "<p>" 
             "$$K_{y}={" + data.KyValue + "}$$"  
             "</p>" 
             
             
             "</div>"
              "</div>" 
               
            ];
           
         
            
        end
        
        function NeDiv=NeSection(~,Nex,Ney,Nez)
            NeDiv=[
           
             "<p class='card-text font-weight-bold'>Condições de Contorno</p>" 
             "<h6>Seção-5.4.3.1-NBR8800-pág:50</h6>"
            
             "<div class='row border border-secondary'>" 
             "<div class='col-md-4 text-center'>"
             "<p class='font-weight-bold'>Ne-x</p>" 
             "<p>"+ Nex.equationString+"</p>" 
            
             "</div>" 
             "<div class='col-md-4 text-center'>" 
             "<p class='font-weight-bold'>Ne-y</p>" 
             "<p>"+ Ney.equationString+"</p>" 
            
             "</div>" 
              "<div class='col-md-4 text-center'>" 
             "<p class='font-weight-bold'>Ne-z</p>" 
               "<p>"+ Nez.equationString+"</p>" 
             "</div>"
              
            
            
              "<div class='col-md-4 text-center'>"
            
             "<p>"+ Nex.equationValue+"</p>" 
            
             "</div>" 
             "<div class='col-md-4 text-center'>" 
             
             "<p>"+ Ney.equationValue+"</p>" 
            
             "</div>" 
              "<div class='col-md-4 text-center'>" 
             
               "<p>"+ Nez.equationValue+"</p>" 
             "</div>"
             
              "<div class='col-md-4 text-center'>"
            
             "<p>Nex="+ Nex.value+"</p>" 
            
             "</div>" 
             "<div class='col-md-4 text-center'>" 
             
             "<p>Ney="+ Ney.value+"</p>" 
            
             "</div>" 
              "<div class='col-md-4 text-center'>" 
             
               "<p>Nez="+ Nez.value+"</p>" 
             "</div>"
             
             
             
             
             "</div>"
             
               
            ];
        end
        
        
        function qaDiv=QaSection(~,esbeltez,esbeltezLim,specificCase)
            if(esbeltez.value <= esbeltezLim.value)
                compare="$$"+esbeltez.value+"<="+esbeltezLim.value+"$$";
                solve=[specificCase.qaValue
                    ];
                
            else
                compare="$$"+esbeltez.value+">="+esbeltezLim.value+"$$";
                solve=[specificCase.befEquationString
                    specificCase.befEquationValue
                    specificCase.befValue
                    specificCase.aefEquationString
                    specificCase.aefValue
                    specificCase.qaEquationString
                    specificCase.qaEquationValue
                    specificCase.qaValue
                    ];
            end
            qaDiv=[
                
            
            "<div class='col-md-12'>"
            "<p class='card-text font-weight-bold'>Fator Qa</p>"
            "<h6>Seção-5.4.3.1-NBR8800-pág:50</h6>"
            "</div>"
            "<div class='row border border-secondary'>"
            
            
            
            "<div class='col-md-6 text-center'>"
            "<p class='font-weight-bold'>Esbeltez</p>"
            "<p>"+ esbeltez.equationString+"</p>"
            
            "</div>"
            "<div class='col-md-6 text-center'>"
            "<p class='font-weight-bold'>Esbeltez Lim</p>"
            "<p>"+ esbeltezLim.equationString+"</p>"
            
            "</div>"
            
            
            
            
            "<div class='col-md-6 text-center'>"
            
            "<p>"+ esbeltez.equationValue+"</p>"
            
            "</div>"
            "<div class='col-md-6 text-center'>"
            
            "<p>"+ esbeltezLim.equationValue+"</p>"
            
            "</div>"
            "<div class='col-md-12 text-center'>"
            "<p>Como:</p>"
            compare
            
            "<p>Logo:</p>"
            solve
            "</div>"
            
            
            
            
            
            
            
            "</div>"
            
            
            ];
        end
        
        
        
         function qsDiv=QsSection(~,esbeltez,esbeltezLim,specificCase)
             
             switch specificCase.compareCase
                 case 1
                     compare=specificCase.compare;
                     solve=[specificCase.qsValue
                         ];
                 case 2
                     compare=specificCase.compare;
                     solve=[
                          specificCase.qsEquationString
                          specificCase.qsEquationValue
                           specificCase.qsValue
                         ];
                     
                 case 3
                     
                      compare=specificCase.compare;
                     solve=[
                          specificCase.qsEquationString
                          specificCase.qsEquationValue
                           specificCase.qsValue
                         ];
                     
                 otherwise
             end
           
            qsDiv=[
                
            
            "<div class='col-md-12'>"
            "<p class='card-text font-weight-bold'>Fator Qs</p>"
            "<h6>Seção-5.4.3.1-NBR8800-pág:50</h6>"
            "</div>"
            "<div class='row border border-secondary'>"
            
            
            
            "<div class='col-md-6 text-center'>"
            "<p class='font-weight-bold'>Esbeltez</p>"
            "<p>"+ esbeltez.equationString+"</p>"
            
            "</div>"
            "<div class='col-md-6 text-center'>"
            "<p class='font-weight-bold'>Esbeltez Lim</p>"
            "<p>"+ esbeltezLim.equationString+"</p>"
            
            "</div>"
            
            
            
            
            "<div class='col-md-6 text-center'>"
            
             "<p>"+ esbeltez.equationValue+"</p>"
            
            "</div>"
            "<div class='col-md-6 text-center'>"
            
%             "<p>"+ esbeltezLim.equationValue+"</p>"
            
            "</div>"
            "<div class='col-md-12 text-center'>"
             "<p>Como:</p>"
            compare
%             
             "<p>Logo:</p>"
             solve
            "</div>"
            
            
            
            
            
            
            
            "</div>"
            
            
            ];
         end
        
        
         function qDiv=QSection(~,Qa,Qs)
             
           
           
            qDiv=[
                
            
            "<div class='col-md-12'>"
            "<p class='card-text font-weight-bold'>Fator Q</p>"
            "<h6>Seção-5.4.3.1-NBR8800-pág:50</h6>"
            "</div>"
            "<div class='row border border-secondary'>"
            
            
            
            "<div class='col-md-12 text-center'>"
            "<p class='font-weight-bold'>Q Resultante</p>"
            "<p>$$Qa*Qs$$</p>"
            
            "</div>"
           
            
            
            
            "<div class='col-md-12 text-center'>"
            
             "<p>$$Qs="+Qa+"*"+Qs+"$$</p>"
            
            "</div>"
           
            "<div class='col-md-12 text-center'>"
            
%             
             "<p>Logo:</p>"
              "<p>$$Qs="+Qa*Qs+"$$</p>"
            "</div>"
            
            
            
            
            
            
            
            "</div>"
            
            
            ];
         end
        
         
         function qsiDiv=QsiSection(~,specificCase,esbeltezReduzida)
             
             switch specificCase.compareCase
                 case 1
                     compare=specificCase.compare;
                     solve=[specificCase.equationString
                            specificCase.value
                         ];
                 case 2
                     compare=specificCase.compare;
                    solve=[specificCase.equationString
                            specificCase.value
                         ];
                     
                
                 otherwise
             end
           
            qsiDiv=[
                
            
            "<div class='col-md-12'>"
            "<p class='card-text font-weight-bold'>Fator X</p>"
            "<h6>Seção-5.4.3.1-NBR8800-pág:50</h6>"
            "</div>"
            "<div class='row border border-secondary'>"
            
            
            
            "<div class='col-md-12 text-center'>"
            "<p class='font-weight-bold'>Esbeltez Reduzida</p>"
            "<p>"+ esbeltezReduzida.equationString+"</p>"
            
            "</div>"
            "<div class='col-md-12 text-center'>"
     
            "<p>"+esbeltezReduzida.equationValue+"</p>"
            
            "</div>"
            
            
            
            
            "<div class='col-md-12 text-center'>"
            
             "<p>$$\gamma_{o}="+esbeltezReduzida.value+"$$</p>"
            
            "</div>"
            "<div class='col-md-12 text-center'>"
            
%             "<p>"+ esbeltezLim.equationValue+"</p>"
            
            "</div>"
            "<div class='col-md-12 text-center'>"
             "<p>Como:</p>"
            compare
%             
             "<p>Logo:</p>"
             solve
            "</div>"
            
            
            
            
            
            
            
            "</div>"
            
            
            ];
         end
         
         
         function nrdDiv=NrdSection(~,nrd)
             
             
           nrdDiv=[
                
            
            "<div class='col-md-12'>"
            "<p class='card-text font-weight-bold'>Nrd</p>"
            "<h6>Seção-5.4.3.1-NBR8800-pág:50</h6>"
            "</div>"
            "<div class='row border border-secondary'>"
            
            
            
            "<div class='col-md-12 text-center'>"
            "<p class='font-weight-bold'>Nrd resultante</p>"
            "<p>"+ nrd.equationString+"</p>"
            
            "</div>"
            "<div class='col-md-12 text-center'>"
     
            "<p>"+nrd.equationValue+"</p>"
            
            "</div>"
            
            
            
            
            "<div class='col-md-12 text-center'>"
            
             "<p>Nrd="+nrd.value+"</p>"
            
            "</div>"
           
            
            
            
            
            
            
            
            "</div>"
            
            
            ];
         end
    end
    
    
    
    
 
end

