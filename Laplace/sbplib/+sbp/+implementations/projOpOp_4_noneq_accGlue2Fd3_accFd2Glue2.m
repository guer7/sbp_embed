function [stencil_g2f,BCU_g2f,HU,M] = projOpOp_4_noneq_accGlue2Fd3_accFd2Glue2
%PROJGLUE_ACC_IO3_BOG2F2_BOF2G1_STENCIL_24_BC_4_24
%    [STENCIL_G2F,BCU_G2F,HU,M] = PROJGLUE_ACC_IO3_BOG2F2_BOF2G1_STENCIL_24_BC_4_24

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    25-Jan-2022 16:34:02

stencil_g2f = [1.59722222222535e-2,1.59722222204315e-2,1.597222220378997e-2,5.179080951254866e-4,-8.958333333343449e-2,-8.263888887995775e-2,-3.958333327803663e-2,-1.399064285159735e-3,5.736111111111812e-1,3.347222222043826e-1,2.361111107424677e-2,1.60765237985226e-3,5.736111111111812e-1,-3.347222222043826e-1,2.361111107424677e-2,-1.60765237985226e-3,-8.958333333343449e-2,8.263888887995775e-2,-3.958333327803663e-2,1.399064285159735e-3,1.59722222222535e-2,-1.59722222204315e-2,1.597222220378997e-2,-5.179080951254866e-4];
if nargout > 1
    BCU_g2f = reshape([1.341740624905692,4.423178224542316e-1,-6.35611494002267e-2,1.729276054378924e-2,-5.53167302855845e-1,1.172159627632877e-1,-3.319137418110962e-3,9.242185093104385e-4,6.38908202848204e-2,-2.387405331218531e-2,1.274676466464043e-2,-2.862207894705228e-3,-2.262441050756662e-3,1.005141899339966e-3,-7.730111566455095e-4,2.866049487240006e-4,-4.331690161033869e-1,6.303811967346724e-1,5.958589701755128e-1,-9.953666763117924e-2,1.563231940958812e-1,-2.65165188672703e-1,2.70163880358656e-1,-6.920773327721712e-2,2.356738707438078e-1,-9.870170997216722e-2,8.326126944361745e-2,-5.529013491588415e-2,-1.078237508240099e-2,4.360386755785679e-3,-2.337691487429867e-3,-1.831461428274805e-4,6.842932952987835e-2,-6.775622850924672e-2,5.091349134491425e-1,5.881755872392744e-1,-1.312136186151238e-1,8.627799277553808e-2,-3.029021467777822e-1,3.371145001571786e-1,1.282975511613846e-1,-5.001767489372975e-2,2.268876375873468e-2,2.353533477722993e-2,-4.431697360217442e-3,1.841464006736656e-3,-1.594736159920649e-3,1.67561248656142e-3,-7.693950530717843e-3,1.123893834685252e-2,-7.186863570716077e-2,5.755595483689876e-1,1.733825339629428e-2,-2.102168242030765e-2,7.978941110015983e-2,-3.40894570638938e-1,-5.325320661673032e-3,1.563730182427575e-2,-3.563128871656478e-2,2.381262199689742e-2,1.127240027705622e-4,-4.780873429763692e-4,1.237438086777438e-3,-1.615787851960477e-3,4.567295834181895e-2,-2.436765427887617e-2,3.863630739249597e-2,-1.013619341898053e-1,1.460864902020493e-3,-8.778134656388555e-5,-1.545418101065409e-2,8.427190475600693e-2,-1.849418124252553e-3,8.934535959535883e-4,1.40467598319978e-2,-3.982484737386233e-2,-1.112956212053364e-4,5.011280370034548e-5,-5.200214850057832e-4,1.433697208738693e-3,-1.497994614328378e-2,8.185925252371372e-3,-8.200405909765104e-3,1.987070566892495e-2,-2.394724022699206e-2,1.196085130544598e-2,-1.084343997385006e-2,-1.163133463908884e-2,-3.306155847613084e-3,1.54960725399915e-3,-1.292253280844808e-3,1.671030793408769e-2,-9.447885040306628e-5,4.257713449454182e-5,-3.349429588401191e-5,-5.125554927718991e-4],[4,24]);
end
if nargout > 2
    HU = [2.1259737557798e-1;1.0260290400758;1.0775123588954;9.8607273802835e-1];
end
if nargout > 3
    M = [2.0;6.666666666666667e-1;4.0e-1;2.857142857142857e-1];
end