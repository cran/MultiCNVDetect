las_bycgd <-
function(Y,X=Y,lambda1=1,lambda2=1){
  eps<-0.0001; alpha<-1;
  lambda1<-max(apply(Y,1,FUN="mad"))*sqrt(dim(Y)[1])
  lambda2<-max(apply(Y,1,FUN="mad"))*sqrt(dim(Y)[1])*sqrt(log(dim(Y)[2]))
  D<-calculate_direction(Y,X,lambda1,lambda2)#����x[k+1]=x[k]+alpha*d�е�d
     
  while( alpha*max(abs(D))>eps ){#��alpha*max��|d|��>epsʱ��BCGD�㷨��������
    alpha<-armijo_rule(Y,X,D,lambda1,lambda2,alpha)#����x[k+1]=x[k]+alpha*d��alpha
    X<-X+alpha*D#x[k+1]=x[k]+alpha*d
    D<-calculate_direction(Y,X,lambda1,lambda2)#������һ�ε����ķ���
   }
  
  return (X)
}