function [stencil_F2C,BC_F2C,HcU,HfU] = intOpOP_orders_2to2_ratio_2to1_accC2F2_accF2C1
%INT_ORDERS_2TO2_RATIO_2TO1_ACCC2F2_ACCF2C1_STENCIL_5_BC_1_3
%    [STENCIL_F2C,BC_F2C,HCU,HFU] = INT_ORDERS_2TO2_RATIO_2TO1_ACCC2F2_ACCF2C1_STENCIL_5_BC_1_3

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    21-May-2018 15:36:08

stencil_F2C = [1.0./4.0,1.0./2.0,1.0./4.0];
if nargout > 1
    BC_F2C = [1.0./2.0,1.0./2.0];
end
if nargout > 2
    HcU = 1.0./2.0;
end
if nargout > 3
    HfU = 1.0./2.0;
end
