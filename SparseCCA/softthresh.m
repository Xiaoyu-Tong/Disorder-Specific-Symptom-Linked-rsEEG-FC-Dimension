function s=softthresh(a,b)
s=sign(a).*max(abs(a)-b,0);