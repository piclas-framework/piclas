function X_new=ChangeBasis(p,X,sVdm_Bezier,dimension)
% ==================================================================
% Description
% ==================================================================
% 

% ==================================================================
% Input
% ==================================================================
%                   ::

% ==================================================================
% Output
% ==================================================================
%                   ::

% ==================================================================
% Function Start
% ==================================================================
switch dimension
  case 1
    %%
    error('this dimension has not yet been implemented in "ChangeBasis(p,P,dimension)"');
  case 2
    %%
    error('this dimension has not yet been implemented in "ChangeBasis(p,P,dimension)"');
  case 3
    %%
    X_temp=zeros(p+1,p+1,3);
    X_new =zeros(p+1,p+1,3);
    for dim=1:3
      for j=1:p+1
        X_temp(:,j,dim)=matmul(sVdm_Bezier,X(:,j,dim));
      end
    end
    
    for dim=1:3
      for i=1:p+1
        X_new(i,:,dim)=matmul(sVdm_Bezier,X_temp(i,:,dim));
      end
    end
    
    
    
    
    
    
  otherwise
    error('wrong dimension in "ChangeBasis(p,P,dimension)"');
end




% ==================================================================
% Function End
% ==================================================================
end