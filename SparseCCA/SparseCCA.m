function [Wx,Wy,R,p]=SparseCCA(X,Y,ncomp,max_iter,lambda1,lambda2)
% % Sparse canonical correlation analysis via penalized matrix decomposition % %
% 
% Input:
% X -- data1: observation x variable
% Y -- data2: observation x variable
% ncomp - number of canonical variates to be learned
% max_iter -- maximal iteration number
% lambda1 -- penality parameter for projection matrix Wx (larger for more sparse)
% lambda2 -- penality parameter for projection matrix Wy (larger for more sparse)
%
% Output:
% Wx -- projection matrix for X
% Wy -- projection matrix for Y
% R -- correlation coefficient
% p -- correlation significance
% 
% by Yu Zhang, Stanford, 2017-10
% yzhangsu@stanford.edu
%
% modified by Xiaoyu Tong, Lehigh, 2022-10
% xit321@lehigh.edu
% 
% Reference: Witten, D.M., Tibshirani, R. and Hastie, T., 2009. 
%               A penalized matrix decomposition, with applications to 
%               sparse principal components and canonical correlation analysis. 
%               Biostatistics, 10(3), pp.515-534.
% 



% standarize
% X=X-repmat(mean(X),size(X,1),1); % no de-mean to ensure the overall
% linearity of transformation for interpretability.
X=X./repmat(std(X),size(X,1),1); % unit variance to ensure sparsity constraint is fair across features
Y=Y-repmat(mean(Y),size(Y,1),1);
Y=Y./repmat(std(Y),size(Y,1),1); % unit variance to ensure sparsity constraint is fair across measures

% sparse cca
Wx=[];
Wy=[];
d=[];
R=[];
p=[];
Xres=X; Yres=Y;
for ncp=1:ncomp
    % initialize
    [u,s,v]=svd(Xres'*Yres);
%     [u,v]=power_svd(Xres, Yres);
    WyInit=v(:,1:ncomp);
    clear u s v;
%     wy=WyInit(:,ncp);
    wy=randn(size(Yres,2),1);
    wyold=randn(length(wy),1);
    wx=randn(size(Xres,2),1);
    for itr=1:max_iter
        if sum(abs(wyold-wy))>1e-6
            % update Wx
            tempwx=Xres'*(Yres*wy);
            lamwx=binarysearch(tempwx,(1-lambda1)*sqrt(size(Xres,2)));
            swx=softthresh(tempwx,lamwx);
            wx=swx./norm(swx);
            
            % update Wy
            wyold=wy;
            tempwy=wx'*Xres'*Yres;
            lamwy=binarysearch(tempwy,(1-lambda2)*sqrt(size(Yres,2)));
            swy=softthresh(tempwy,lamwy);
            wy=swy'./norm(swy);
        else
            break
        end
    end
    
    Wx(:,ncp)=wx;
    Wy(:,ncp)=wy;
    d(ncp)=sum((Xres*wx)'*(Yres*wy));
    
    Xres=[Xres;sqrt(d(ncp))*wx'];
    Yres=[Yres;-sqrt(d(ncp))*wy'];
    
    [R(ncp),p(ncp)]=corr(X*wx,Y*wy);
%     R(ncp)=(wx'*X'*Y*wy)/sqrt((wx'*X'*X*wx)*(wy'*Y'*Y*wy));
    a =0;
end



