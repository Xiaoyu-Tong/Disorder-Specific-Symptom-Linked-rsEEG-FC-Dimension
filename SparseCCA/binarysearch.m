function lam=binarysearch(arg,sumabs)

if norm(arg,2) == 0 || (sum(abs(arg/norm(arg, 2))) - sumabs)<1e-6
    lam = 0;
    return;
end
lam1=0;
lam2=max(abs(arg))-1e-5;
iter=1;
while(iter <= 150)
    sh=softthresh(arg,(lam1+lam2)/2);
    if(sum(abs(sh/norm(sh)))<sumabs)
      lam2=(lam1+lam2)/2;
    else
      lam1=(lam1+lam2)/2;
    end
    if((lam2-lam1)<1e-6) 
        break
    end
    iter=iter+1;
end
lam=(lam1+lam2)/2;
if iter==150
    warning('Did not quite converge');
end




