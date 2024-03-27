function F=Wfield(X,vec)
global eps Rm bp bm gam
 v1=vec(1);
 v2=vec(2);
 F(1)=v2;
 F(2)=Wp(v1);
 if norm(F)>10
     F=F*exp(10-norm(F));
 end
 F=transpose(F);
end

