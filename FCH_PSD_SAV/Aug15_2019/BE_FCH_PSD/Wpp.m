function Wpp=Wpp(u)
  global bm bp gam
  Wpp=(2*u-bm-bp).*(u-(bp+bm-gam)/2)+(u-bm).*(u-bp); 
end

