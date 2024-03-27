function [Dv] = CalcDerivative(v, delta)
for ind=1:1:length(v)
    if ind==1
        Dv(ind) = (v(ind+1)-v(ind))/delta;
    elseif ind == length(v)
        Dv(ind) = (v(ind)-v(ind-1))/delta;
    else
        Dv(ind) = (v(ind+1)-v(ind-1))/(2*delta);
    end
end
end
