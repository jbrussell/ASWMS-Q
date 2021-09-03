function [AngDiff] = angdiff(a,b)
%difference between two angles (a - b) in radians
AngDiff = [];
for ii = 1:length(a)
    AngDiff(ii,:) = diff(unwrap([b(ii),a(ii)]));
end
end

