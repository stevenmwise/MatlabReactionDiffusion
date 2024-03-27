function F=LGfield(t,vec)
global eps eta a bm bp Xbl Phibl lamb
 v1=vec(1);
 v2=vec(2);
 if (abs(t)<=Xbl(end))
    YI=interp1(Xbl,Phibl,abs(t));  %Phibl is even, defined for Xbl>0, but t<0
 else
    YI=b;
 end
 F(1)=v2;
 F(2)=(Wpp(YI)+lamb)*v1; 
 if norm(F)>10
     F=F*exp(10-norm(F));
 end
 F=transpose(F);
end
