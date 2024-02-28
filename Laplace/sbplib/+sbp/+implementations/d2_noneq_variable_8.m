function [H, HI, D1, D2, DI] = d2_noneq_variable_8(N, h, options)
    % N: Number of grid points
    % h: grid spacing
    % options: struct containing options for constructing the operator
    %          current options are: 
    %               options.stencil_type ('minimal','nonminimal','wide')
    %               options.AD ('upwind', 'op')

    % BP: Number of boundary points
    % order: Accuracy of interior stencil
    BP = 8;
    order = 8;
    if(N<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    %%%% Norm matrix %%%%%%%%
    P = zeros(BP, 1);
    P0 = 1.0758368078310e-01;
    P1 = 6.1909685107891e-01;
    P2 = 9.6971176519117e-01;
    P3 = 1.1023441350947e+00;
    P4 = 1.0244688965833e+00;
    P5 = 9.9533550116831e-01;
    P6 = 1.0008236941028e+00;
    P7 = 9.9992060631812e-01;

    for i = 0:BP - 1
        P(i + 1) = eval(['P' num2str(i)]);
    end

    Hv = ones(N, 1);
    Hv(1:BP) = P;
    Hv(end - BP + 1:end) = flip(P);
    Hv = h * Hv;
    H = spdiags(Hv, 0, N, N);
    HI = spdiags(1 ./ Hv, 0, N, N);
    %%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Q matrix %%%%%%%%%%%

    % interior stencil
    d = [1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280];
    d = repmat(d, N, 1);
    Q = spdiags(d, -order / 2:order / 2, N, N);

    % Boundaries
    Q0_0 = -5.0000000000000e-01;
    Q0_1 = 6.7284756079369e-01;
    Q0_2 = -2.5969732837062e-01;
    Q0_3 = 1.3519390385721e-01;
    Q0_4 = -6.9678474730984e-02;
    Q0_5 = 2.6434024071371e-02;
    Q0_6 = -5.5992311465618e-03;
    Q0_7 = 4.9954552590464e-04;
    Q0_8 = 0.0000000000000e+00;
    Q0_9 = 0.0000000000000e+00;
    Q0_10 = 0.0000000000000e+00;
    Q0_11 = 0.0000000000000e+00;
    Q1_0 = -6.7284756079369e-01;
    Q1_1 = 0.0000000000000e+00;
    Q1_2 = 9.4074021172233e-01;
    Q1_3 = -4.0511642426516e-01;
    Q1_4 = 1.9369192209331e-01;
    Q1_5 = -6.8638079843479e-02;
    Q1_6 = 1.3146457241484e-02;
    Q1_7 = -9.7652615479254e-04;
    Q1_8 = 0.0000000000000e+00;
    Q1_9 = 0.0000000000000e+00;
    Q1_10 = 0.0000000000000e+00;
    Q1_11 = 0.0000000000000e+00;
    Q2_0 = 2.5969732837062e-01;
    Q2_1 = -9.4074021172233e-01;
    Q2_2 = 0.0000000000000e+00;
    Q2_3 = 9.4316393361096e-01;
    Q2_4 = -3.5728039257451e-01;
    Q2_5 = 1.1266686855013e-01;
    Q2_6 = -1.8334941452280e-02;
    Q2_7 = 8.2741521740941e-04;
    Q2_8 = 0.0000000000000e+00;
    Q2_9 = 0.0000000000000e+00;
    Q2_10 = 0.0000000000000e+00;
    Q2_11 = 0.0000000000000e+00;
    Q3_0 = -1.3519390385721e-01;
    Q3_1 = 4.0511642426516e-01;
    Q3_2 = -9.4316393361096e-01;
    Q3_3 = 0.0000000000000e+00;
    Q3_4 = 8.7694387866575e-01;
    Q3_5 = -2.4698058719506e-01;
    Q3_6 = 4.7291642094198e-02;
    Q3_7 = -4.0135203618880e-03;
    Q3_8 = 0.0000000000000e+00;
    Q3_9 = 0.0000000000000e+00;
    Q3_10 = 0.0000000000000e+00;
    Q3_11 = 0.0000000000000e+00;
    Q4_0 = 6.9678474730984e-02;
    Q4_1 = -1.9369192209331e-01;
    Q4_2 = 3.5728039257451e-01;
    Q4_3 = -8.7694387866575e-01;
    Q4_4 = 0.0000000000000e+00;
    Q4_5 = 8.1123946853807e-01;
    Q4_6 = -2.0267150541446e-01;
    Q4_7 = 3.8680398901392e-02;
    Q4_8 = -3.5714285714286e-03;
    Q4_9 = 0.0000000000000e+00;
    Q4_10 = 0.0000000000000e+00;
    Q4_11 = 0.0000000000000e+00;
    Q5_0 = -2.6434024071371e-02;
    Q5_1 = 6.8638079843479e-02;
    Q5_2 = -1.1266686855013e-01;
    Q5_3 = 2.4698058719506e-01;
    Q5_4 = -8.1123946853807e-01;
    Q5_5 = 0.0000000000000e+00;
    Q5_6 = 8.0108544742793e-01;
    Q5_7 = -2.0088756283071e-01;
    Q5_8 = 3.8095238095238e-02;
    Q5_9 = -3.5714285714286e-03;
    Q5_10 = 0.0000000000000e+00;
    Q5_11 = 0.0000000000000e+00;
    Q6_0 = 5.5992311465618e-03;
    Q6_1 = -1.3146457241484e-02;
    Q6_2 = 1.8334941452280e-02;
    Q6_3 = -4.7291642094198e-02;
    Q6_4 = 2.0267150541446e-01;
    Q6_5 = -8.0108544742793e-01;
    Q6_6 = 0.0000000000000e+00;
    Q6_7 = 8.0039405922650e-01;
    Q6_8 = -2.0000000000000e-01;
    Q6_9 = 3.8095238095238e-02;
    Q6_10 = -3.5714285714286e-03;
    Q6_11 = 0.0000000000000e+00;
    Q7_0 = -4.9954552590464e-04;
    Q7_1 = 9.7652615479254e-04;
    Q7_2 = -8.2741521740941e-04;
    Q7_3 = 4.0135203618880e-03;
    Q7_4 = -3.8680398901392e-02;
    Q7_5 = 2.0088756283071e-01;
    Q7_6 = -8.0039405922650e-01;
    Q7_7 = 0.0000000000000e+00;
    Q7_8 = 8.0000000000000e-01;
    Q7_9 = -2.0000000000000e-01;
    Q7_10 = 3.8095238095238e-02;
    Q7_11 = -3.5714285714286e-03;

    for i = 1:BP

        for j = 1:BP
            Q(i, j) = eval(['Q' num2str(i - 1) '_' num2str(j - 1)]);
            Q(N + 1 - i, N + 1 - j) = -eval(['Q' num2str(i - 1) '_' num2str(j - 1)]);
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Undivided difference operators %%%%
    % Closed with zeros at the first boundary nodes.
    m = N;

    DD_4 = (diag(ones(m - 2, 1), 2) - 4 * diag(ones(m - 1, 1), 1) + 6 * diag(ones(m, 1), 0) - 4 * diag(ones(m - 1, 1), -1) + diag(ones(m - 2, 1), -2));
    DD_4(1:6, 1:8) = [0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0.70921010190504348684e1 -0.14196080536361841322e2 0.11072881931325435634e2 -0.50473576941871051066e1 0.10784552801730759259e1 0 0 0; 0 0.13740993382151221352e1 -0.42105600869792757010e1 0.54761010136211975317e1 -0.35797005751940657417e1 0.94006031033702177578e0 0 0; 0 0 0.82467928104463767301e0 -0.33274694995849432461e1 0.52587584638857303123e1 -0.37020511582893568152e1 0.94608291294393207601e0 0; 0 0 0 0.86436129166612654748e0 -0.37325441295306179390e1 0.57924699560798105338e1 -0.39066885960487908497e1 0.98240147783347170744e0; ];
    DD_4(m - 5:m, m - 7:m) = [0.98240147783347170744e0 -0.39066885960487908497e1 0.57924699560798105338e1 -0.37325441295306179390e1 0.86436129166612654748e0 0 0 0; 0 0.94608291294393207601e0 -0.37020511582893568152e1 0.52587584638857303123e1 -0.33274694995849432461e1 0.82467928104463767301e0 0 0; 0 0 0.94006031033702177578e0 -0.35797005751940657417e1 0.54761010136211975317e1 -0.42105600869792757010e1 0.13740993382151221352e1 0; 0 0 0 0.10784552801730759259e1 -0.50473576941871051066e1 0.11072881931325435634e2 -0.14196080536361841322e2 0.70921010190504348684e1; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; ];
    DD_4 = sparse(DD_4);

    DD_5 = (-diag(ones(m - 3, 1), -3) + 5 * diag(ones(m - 2, 1), -2) - 10 * diag(ones(m - 1, 1), -1) + 10 * diag(ones(m, 1), 0) - 5 * diag(ones(m - 1, 1), 1) + diag(ones(m - 2, 1), 2));
    DD_5(1:7, 1:9) = [0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; -0.82098088052411907132e1 0.18024024120655590699e2 -0.17692096674769448317e2 0.12181944917152947702e2 -0.53922764008653796295e1 0.10882128430674802600e1 0 0 0; 0 -0.13913240333052950751e1 0.50983574122643719988e1 -0.89139255753021486176e1 0.89492514379851643542e1 -0.47003015516851088789e1 0.95794231004301621873e0 0 0; 0 0 -0.80388595981635823711e0 0.40861386922975635215e1 -0.87645974398095505205e1 0.92551278957233920380e1 -0.47304145647196603800e1 0.95763137632461357808e0 0; 0 0 0 -0.85214912336218661144e0 0.46656801619132724238e1 -0.96541165934663508896e1 0.97667214901219771242e1 -0.49120073891673585372e1 0.98587145396064649030e0; ];
    DD_5(m - 5:m, m - 8:m) = [-0.98587145396064649030e0 0.49120073891673585372e1 -0.97667214901219771242e1 0.96541165934663508896e1 -0.46656801619132724238e1 0.85214912336218661144e0 0 0 0; 0 -0.95763137632461357808e0 0.47304145647196603800e1 -0.92551278957233920380e1 0.87645974398095505205e1 -0.40861386922975635215e1 0.80388595981635823711e0 0 0; 0 0 -0.95794231004301621873e0 0.47003015516851088789e1 -0.89492514379851643542e1 0.89139255753021486176e1 -0.50983574122643719988e1 0.13913240333052950751e1 0; 0 0 0 -0.10882128430674802600e1 0.53922764008653796295e1 -0.12181944917152947702e2 0.17692096674769448317e2 -0.18024024120655590699e2 0.82098088052411907132e1; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; ];
    DD_5 = sparse(DD_5);

    DD_6 = (diag(ones(m - 3, 1), 3) - 6 * diag(ones(m - 2, 1), 2) + 15 * diag(ones(m - 1, 1), 1) - 20 * diag(ones(m, 1), 0) + 15 * diag(ones(m - 1, 1), -1) - 6 * diag(ones(m - 2, 1), -2) + diag(ones(m - 3, 1), -3));
    DD_6(1:7, 1:10) = [0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0; 0.92604272237009487731e1 -0.21899951980342041757e2 0.25706973995952182034e2 -0.23795532642767976509e2 0.16176829202596138889e2 -0.65292770584048815597e1 0.10805312592656301300e1 0 0 0; 0 0.14058275749850312569e1 -0.59637699688344396810e1 0.13135580487711615439e2 -0.17898502875970328708e2 0.14100904655055326637e2 -0.57476538602580973124e1 0.96761398731089236944e0 0 0; 0 0 0.78692381135906040550e0 -0.48340889923923043812e1 0.13146896159714325781e2 -0.18510255791446784076e2 0.14191243694158981140e2 -0.57457882579476814685e1 0.96506937655440259938e0 0; 0 0 0 0.84209241882516286974e0 -0.55988161942959269085e1 0.14481174890199526334e2 -0.19533442980243954248e2 0.14736022167502075612e2 -0.59152287237638789418e1 0.98819842177699528293e0; ];
    DD_6(m - 6:m, m - 9:m) = [0.98819842177699528293e0 -0.59152287237638789418e1 0.14736022167502075612e2 -0.19533442980243954248e2 0.14481174890199526334e2 -0.55988161942959269085e1 0.84209241882516286974e0 0 0 0; 0 0.96506937655440259938e0 -0.57457882579476814685e1 0.14191243694158981140e2 -0.18510255791446784076e2 0.13146896159714325781e2 -0.48340889923923043812e1 0.78692381135906040550e0 0 0; 0 0 0.96761398731089236944e0 -0.57476538602580973124e1 0.14100904655055326637e2 -0.17898502875970328708e2 0.13135580487711615439e2 -0.59637699688344396810e1 0.14058275749850312569e1 0; 0 0 0 0.10805312592656301300e1 -0.65292770584048815597e1 0.16176829202596138889e2 -0.23795532642767976509e2 0.25706973995952182034e2 -0.21899951980342041757e2 0.92604272237009487731e1; 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0; ];
    DD_6 = sparse(DD_6);

    DD_7 = (-diag(ones(m - 4, 1), -4) + 7 * diag(ones(m - 3, 1), -3) - 21 * diag(ones(m - 2, 1), -2) + 35 * diag(ones(m - 1, 1), -1) - 35 * diag(ones(m, 1), 0) + 21 * diag(ones(m - 1, 1), 1) - 7 * diag(ones(m - 2, 1), 2) + diag(ones(m - 3, 1), 3));
    DD_7(1:8, 1:11) = [0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0; -0.10257962606384161889e2 0.25816283570514673553e2 -0.35082323899232605336e2 0.40909341259657284245e2 -0.37745934806057657407e2 0.22852469704417085459e2 -0.75637188148594109098e1 0.10718455919447922848e1 0 0 0; 0 -0.14183700945143243060e1 0.68109221539151178574e1 -0.18129991367652287552e2 0.31322380032948075240e2 -0.32902110861795762152e2 0.20116788510903340593e2 -0.67732979111762465861e1 0.97367953737208690559e0 0 0; 0 0 -0.77264857229409862775e0 0.55732122985135573041e1 -0.18405654623600056093e2 0.32392947635031872133e2 -0.33112901953037622660e2 0.20110258902816885140e2 -0.67554856358808181957e1 0.97027194845028099974e0 0; 0 0 0 -0.83355973075426177031e0 0.65319522266785813933e1 -0.20273644846279336868e2 0.34183525215426919935e2 -0.34384051724171509761e2 0.20703300533173576296e2 -0.69173889524389669805e1 0.98986727836499775519e0; ];
    DD_7(m - 6:m, m - 10:m) = [-0.98986727836499775519e0 0.69173889524389669805e1 -0.20703300533173576296e2 0.34384051724171509761e2 -0.34183525215426919935e2 0.20273644846279336868e2 -0.65319522266785813933e1 0.83355973075426177031e0 0 0 0; 0 -0.97027194845028099974e0 0.67554856358808181957e1 -0.20110258902816885140e2 0.33112901953037622660e2 -0.32392947635031872133e2 0.18405654623600056093e2 -0.55732122985135573041e1 0.77264857229409862775e0 0 0; 0 0 -0.97367953737208690559e0 0.67732979111762465861e1 -0.20116788510903340593e2 0.32902110861795762152e2 -0.31322380032948075240e2 0.18129991367652287552e2 -0.68109221539151178574e1 0.14183700945143243060e1 0; 0 0 0 -0.10718455919447922848e1 0.75637188148594109098e1 -0.22852469704417085459e2 0.37745934806057657407e2 -0.40909341259657284245e2 0.35082323899232605336e2 -0.25816283570514673553e2 0.10257962606384161889e2; 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0; ];
    DD_7 = sparse(DD_7);

    DD_8 = (diag(ones(m - 4, 1), 4) - 8 * diag(ones(m - 3, 1), 3) + 28 * diag(ones(m - 2, 1), 2) - 56 * diag(ones(m - 1, 1), 1) + 70 * diag(ones(m, 1), 0) - 56 * diag(ones(m - 1, 1), -1) + 28 * diag(ones(m - 2, 1), -2) - 8 * diag(ones(m - 3, 1), -3) + diag(ones(m - 4, 1), -4));
    DD_8(1:8, 1:12) = [0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0.11211983054345146839e2 -0.29767555907566203748e2 0.45789440151310044599e2 -0.64530162797168947524e2 0.75491869612115314813e2 -0.60939919211778894557e2 0.30254875259437643639e2 -0.85747647355583382784e1 0.10642345748642342173e1 0 0 0; 0 0.14294303785648149613e1 -0.76427065234692260122e1 0.23888038475126051744e2 -0.50115808052716920383e2 0.65804221723591524304e2 -0.53644769362408908249e2 0.27093191644704986344e2 -0.77894362989766952447e1 0.97783801558437253538e0 0 0; 0 0 0.76035645561041265933e0 -0.63048462739199410204e1 0.24540872831466741457e2 -0.51828716216050995413e2 0.66225803906075245320e2 -0.53627357074178360373e2 0.27021942543523272783e2 -0.77621755876022479980e1 0.97411941507587258409e0 0; 0 0 0 0.82615990808320718845e0 -0.74650882590612358780e1 0.27031526461705782491e2 -0.54693640344683071895e2 0.68768103448343019521e2 -0.55208801421796203457e2 0.27669555809755867922e2 -0.79189382269199820415e1 0.99112262457261614959e0; ];
    DD_8(m - 7:m, m - 11:m) = [0.99112262457261614959e0 -0.79189382269199820415e1 0.27669555809755867922e2 -0.55208801421796203457e2 0.68768103448343019521e2 -0.54693640344683071895e2 0.27031526461705782491e2 -0.74650882590612358780e1 0.82615990808320718845e0 0 0 0; 0 0.97411941507587258409e0 -0.77621755876022479980e1 0.27021942543523272783e2 -0.53627357074178360373e2 0.66225803906075245320e2 -0.51828716216050995413e2 0.24540872831466741457e2 -0.63048462739199410204e1 0.76035645561041265933e0 0 0; 0 0 0.97783801558437253538e0 -0.77894362989766952447e1 0.27093191644704986344e2 -0.53644769362408908249e2 0.65804221723591524304e2 -0.50115808052716920383e2 0.23888038475126051744e2 -0.76427065234692260122e1 0.14294303785648149613e1 0; 0 0 0 0.10642345748642342173e1 -0.85747647355583382784e1 0.30254875259437643639e2 -0.60939919211778894557e2 0.75491869612115314813e2 -0.64530162797168947524e2 0.45789440151310044599e2 -0.29767555907566203748e2 0.11211983054345146839e2; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; ];
    DD_8 = sparse(DD_8);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%% Difference operators %%
    D1 = H \ Q;

    % Helper functions for constructing D2(c)
    % TODO: Consider changing sparse(diag(...)) to spdiags(....)
    min_inds =  sbp.implementations.d2_sparsity_pattern_inds(m,order,BP,0,4);
    nonmin_inds = sbp.implementations.d2_sparsity_pattern_inds(m,order,BP,1,4);
    % Minimal 9 point stencil width
    function D2 = D2_fun_minimal(c)
        % Here we add variable diffusion
        C1 = sparse(diag(c));
        C2 = 1/2 * diag(ones(m - 1, 1), -1) + 1/2 * diag(ones(m, 1), 0); C2(1, 2) = 1/2;
        C3 = 3/10 * diag(ones(m - 1, 1), -1) + 3/10 * diag(ones(m - 1, 1), 1) + 2/5 * diag(ones(m, 1), 0); C3(1, 3) = 3/10; C3(m, m - 2) = 3/10;
        C4 = 1/4 * diag(ones(m - 2, 1), -2) + 1/4 * diag(ones(m - 1, 1), -1) + 1/4 * diag(ones(m - 1, 1), 1) + 1/4 * diag(ones(m, 1), 0); C4(2, 4) = 1/4; C4(1, 3) = 1/4; C4(1, 4) = 1/4; C4(m, m - 2) = 1/4; C4(m, m - 3) = 1/4;

        C2 = sparse(diag(C2 * c));
        C3 = sparse(diag(C3 * c));
        C4 = sparse(diag(C4 * c));

        % Remainder term added to wide second derivative operator
        R = (1/78400 / h) * transpose(DD_8) * C1 * DD_8 + (1/14700 / h) * transpose(DD_7) * C2 * DD_7 + (1/2520 / h) * transpose(DD_6) * C3 * DD_6 + (1/350 / h) * transpose(DD_5) * C4 * DD_5;

        D2 = D1 * C1 * D1 - H \ R;

        % Remove potential round off zeros
        D2_tmp = sparse(m,m);
        D2_tmp(min_inds) = D2(min_inds);
        D2 = D2_tmp;
    end

    %  Non-minimal 11 point stencil width
    function D2 = D2_fun_nonminimal(c)
        % Here we add variable diffusion
        C1 = sparse(diag(c));
        C2 = 1/2 * diag(ones(m - 1, 1), -1) + 1/2 * diag(ones(m, 1), 0); C2(1, 2) = 1/2;
        C3 = 3/10 * diag(ones(m - 1, 1), -1) + 3/10 * diag(ones(m - 1, 1), 1) + 2/5 * diag(ones(m, 1), 0); C3(1, 3) = 3/10; C3(m, m - 2) = 3/10;

        C2 = sparse(diag(C2 * c));
        C3 = sparse(diag(C3 * c));

        % Remainder term added to wide second derivative operator
        R = (1/78400 / h) * transpose(DD_8) * C1 * DD_8 + (1/14700 / h) * transpose(DD_7) * C2 * DD_7 + (1/2520 / h) * transpose(DD_6) * C3 * DD_6;
        D2 = D1 * C1 * D1 - H \ R;
        
        % Remove potential round off zeros
        D2_tmp = sparse(m,m);
        D2_tmp(nonmin_inds) = D2(nonmin_inds);
        D2 = D2_tmp;
    end

    % Wide stencil
    function D2 = D2_fun_wide(c)
        % Here we add variable diffusion
        C1 = sparse(diag(c));
        D2 = D1 * C1 * D1;
    end

    switch options.stencil_width
        case 'minimal'
            D2 = @D2_fun_minimal;
        case 'nonminimal'
            D2 = @D2_fun_nonminimal;
        case 'wide'
            D2 = @D2_fun_wide;
        otherwise
            error('No option %s for stencil width', options.stencil_width)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Artificial dissipation operator %%%
    switch options.AD
        case 'upwind'
            % This is the choice that yield 7th order Upwind
            DI = H \ (transpose(DD_4) * DD_4) * (-1/280);
        case 'op'
            % This choice will preserve the order of the underlying
            % Non-dissipative D1 SBP operator
            DI = H \ (transpose(DD_5) * DD_5) * (-1 / (5 * 280));
            % Notice that you can use any negative number instead of (-1/(5*280))
        otherwise
            error("Artificial dissipation options '%s' not implemented.", option.AD)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
end