function initrange = psocheckpopulationinitrange(initrange,LB,UB)
% Automatically adjust PopInitRange according to provided LB and UB.

% if size(LB,2) == 1
%     if ~isinf(LB), lowerRange = LB; end
%     if ~isinf(UB), upperRange = UB; end
% else
    lowerRange = initrange(1,:) ;
    upperRange = initrange(2,:) ;
    
    lowerInf = isinf(LB) ;
    index = false(size(lowerRange,2),1) ;
    index(~lowerInf) = LB(~lowerInf) ~= lowerRange(~lowerInf) ;
    lowerRange(index) = LB(index) ;
    
    upperInf = isinf(UB) ;
    index = false(size(upperRange,2),1) ;
    index(~upperInf) = UB(~upperInf) ~= upperRange(~upperInf) ;
    upperRange(index) = UB(index) ;
% end

initrange = [lowerRange; upperRange] ;