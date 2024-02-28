function [stencil_F2C,BC_F2C,HcU,HfU] = intOpOP_orders_4to4_ratio_2to1_accC2F2_accF2C3
%INT_ORDERS_4TO4_RATIO_2TO1_ACCC2F2_ACCF2C3_STENCIL_9_BC_3_11
%    [STENCIL_F2C,BC_F2C,HCU,HFU] = INT_ORDERS_4TO4_RATIO_2TO1_ACCC2F2_ACCF2C3_STENCIL_9_BC_3_11

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    21-May-2018 15:36:01

stencil_F2C = [7.0./2.56e2,-1.0./3.2e1,-7.0./6.4e1,9.0./3.2e1,8.5e1./1.28e2,9.0./3.2e1,-7.0./6.4e1,-1.0./3.2e1,7.0./2.56e2];
if nargout > 1
    BC_F2C = reshape([7.523257802630956e-2,2.447812262221267e-1,-1.679313063616916e-1,1.290510950666589,6.315723344677289e-3,1.67178747937954e-1,1.982667903557025,-7.554383893379468e-1,7.215271362899867e-1,-2.820478831383137,1.807034411305536,-7.589683751979843e-1,-7.685973268458095e-1,3.965751544535173e-1,4.119789638051451e-1,1.574000556785898e-1,-1.113618964927466e-1,3.631000220124657e-1,1.639694219982991,-9.114074742274862e-1,4.952715520862987e-1,-4.524162151353456e-1,2.65481378213589e-1,-2.268273408674622e-1,3.455365956773523e-1,-2.009765975089827e-1,1.722179564305811e-1,-5.52933488235368e-1,3.18109976271229e-1,-2.171481232558432e-1,1.033835580108027e-1,-5.911351224351343e-2,3.960076712054991e-2],[3,11]);
end
if nargout > 2
    t2 = [1.7e1./4.8e1,5.9e1./4.8e1,4.3e1./4.8e1,4.9e1./4.8e1];
    HcU = t2;
end
if nargout > 3
    HfU = t2;
end
