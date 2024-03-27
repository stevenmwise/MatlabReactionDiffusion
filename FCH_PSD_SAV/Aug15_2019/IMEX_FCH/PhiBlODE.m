  function F=PhiBlODE(t,x)
  global Xbl bm bp Phibl
        v1=x(1);
        v2=x(2);
        if (abs(t)<=Xbl(end))
            YI=interp1(Xbl,Phibl,abs(t));  %Phibl is even, defined for Xbl>0, but t<0
        end
        F(1)=v2;
        F(2)=Wpp(YI)*v1+((Wpp(bm)-Wpp(YI))/Wpp(bm));
        F=transpose(F);
    end
