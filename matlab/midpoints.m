function [xmid, ymid] = midpoints(x,y)
xmid = diag((x + x')/2,1);
ymid = diag((y + y')/2,1);
end
