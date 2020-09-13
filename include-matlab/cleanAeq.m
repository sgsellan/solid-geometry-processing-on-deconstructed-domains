function [Aeq,b1] = cleanAeq(Aeq,csgvec)
    remove = [];
    b1 = [];
    for i=1:size(Aeq,1)
       vec = abs(Aeq(i,:))';
       nonzero = find(vec>0);
       if any(ismember(nonzero,find(csgvec==.5)))
           outside = find((vec>0).*(csgvec==0));
           b1 = [b1;outside];
           if any(outside)
               remove = [remove,i];
           end
       end
%         ss = sum(ismember(nonzero,find(csgvec==1)));
%         if ss<3
%             remove = [remove,i];
%         end
    end
    Aeq(remove,:) = [];
end