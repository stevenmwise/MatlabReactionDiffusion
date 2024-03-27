function Wp=Wp(u)
  global bm bp gam
Wp=(u-bm).*(u-bp).*(u-(bp+bm-gam)/2);
end

