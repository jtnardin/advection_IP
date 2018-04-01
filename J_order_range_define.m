function range = J_order_range_define(num_meth,IC,sigma,N)

    if strcmp(IC,'_gauss')
        
        if strcmp(num_meth,'upwind')

          if sigma <= 2
              range = 1:5;
          elseif sigma <= 4
              range = 1:4;
          else
              range = 1:7;
          end

        elseif strcmp(num_meth,'laxwend')


          if sigma <= 2
              range = 1:4;
          elseif sigma <= 4
              range = 1:3;
          else
              range = 1:7;
          end

        elseif strcmp(num_meth,'beamwarm')

          if sigma <= 3
              range = 1:4;
          elseif sigma <= 4
              range = 1:3;
          else
              range = 1:7;
          end

        elseif strcmp(num_meth,'upwindfl')
            
          if sigma <= 4
              range = 1:4;
          
          else
              range = 1:7;
          end
            
          
        end
    
    
        
        
    elseif strcmp(IC,'_front')
        
        
        
        if strcmp(num_meth,'upwind')

               if sigma <= 4
                   range = 1:7;
               else
                   range = 1:5;
               end

           
           
        elseif strcmp(num_meth,'laxwend')
           
           if N == 1
               if sigma <= 3
                   range = 2:6;
               else
                   range = 1:4;
               end
               
           else
               
               if sigma <= 4
                   range = 1:7;
               else
                   range = 1:4;
               end
               
           end
           
           
       elseif strcmp(num_meth,'beamwarm')
          
           if sigma <= 2
               range = 1:5;
           elseif sigma <= 4
               range = 1:4;
           elseif sigma <= 6
               range = 1:3;
           else
               range = 1:7;
           end
              
           
       elseif strcmp(num_meth,'upwindfl')
          
           if sigma <= 1
               range = 1:7;
           elseif sigma <= 4
               range = 1:4;
           elseif sigma <= 6
               range = 1:4;
           else
               range = 1:7;
           end
           
       end
        
        
        
        
    end



end