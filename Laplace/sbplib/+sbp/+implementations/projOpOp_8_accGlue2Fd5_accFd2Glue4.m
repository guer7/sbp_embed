function [stencil_g2f,BCU_g2f,HU,M] = projOpOp_8_accGlue2Fd5_accFd2Glue4
%PROJGLUE_TRAD_IO7_BOG2F4_BOF2G3_STENCIL_80_BC_8_96
%    [STENCIL_G2F,BCU_G2F,HU,M] = PROJGLUE_TRAD_IO7_BOG2F4_BOF2G3_STENCIL_80_BC_8_96

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    15-Nov-2021 12:21:27

mt1 = [7.408781362638712e-4,-1.142669802121755e-3,-1.802806353580475e-4,-6.37404966582161e-4,-1.762961779466425e-3,6.045002015071175e-5];
mt2 = [-2.223073436426046e-5,3.697268850065437e-5,-6.765181342333469e-3,1.024048776977475e-2,1.594305718053289e-3,5.756235428802924e-3];
mt3 = [1.23314631758229e-2,-5.446475748149276e-4,1.556729530517112e-4,-3.32750347465172e-4,3.035162355661841e-2,-4.053371068669698e-2];
mt4 = [-6.85462064724485e-3,-2.321598903583071e-2,-3.517982040402527e-2,2.183388732586487e-3,-4.449037497908333e-4,1.330989842874565e-3];
mt5 = [-9.953923064518228e-2,9.080404094268676e-2,2.113416732033176e-2,5.580527436384806e-2,4.918029995180156e-2,-5.105378220044177e-3];
mt6 = [6.229808747045193e-4,-3.105625004885142e-3,5.752119102946365e-1,3.552371110275669e-2,-1.569357175578241e-2,-8.589288121477236e-2];
mt7 = [-2.456898094413272e-2,7.667673839548403e-3,-3.115193436010823e-4,4.658424036032496e-3,5.752119102946365e-1,-3.552371110275669e-2];
mt8 = [-1.569357175578241e-2,8.589288121477236e-2,-2.456898094413272e-2,-7.667673839548403e-3,-3.115193436010823e-4,-4.658424036032496e-3];
mt9 = [-9.953923064518228e-2,-9.080404094268676e-2,2.113416732033176e-2,-5.580527436384806e-2,4.918029995180156e-2,5.105378220044177e-3];
mt10 = [6.229808747045193e-4,3.105625004885142e-3,3.035162355661841e-2,4.053371068669698e-2,-6.85462064724485e-3,2.321598903583071e-2];
mt11 = [-3.517982040402527e-2,-2.183388732586487e-3,-4.449037497908333e-4,-1.330989842874565e-3,-6.765181342333469e-3,-1.024048776977475e-2];
mt12 = [1.594305718053289e-3,-5.756235428802924e-3,1.23314631758229e-2,5.446475748149276e-4,1.556729530517112e-4,3.32750347465172e-4];
mt13 = [7.408781362638712e-4,1.142669802121755e-3,-1.802806353580475e-4,6.37404966582161e-4,-1.762961779466425e-3,-6.045002015071175e-5];
mt14 = [-2.223073436426046e-5,-3.697268850065437e-5];
stencil_g2f = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9,mt10,mt11,mt12,mt13,mt14],1,80);
if nargout > 1
    mt15 = [1.604311924593926,3.068012904657915e-1,5.237010798122134e-1,-2.70179385209028e-2,-7.883209999863861e-2,8.212954060596994e-4];
    mt16 = [-4.92505400397121e-3,8.515716628527629e-3,-1.429967156364542e-1,-9.478306783936974e-2,4.202003135646262e-1,1.177188806529701e-1];
    mt17 = [-1.395276459035892e-1,-9.728177704174801e-2,3.258917921469439e-2,1.86213532534563e-2,4.349488749493367e-2,-8.88670865257922e-3];
    mt18 = [-3.525329836822257e-2,-3.874736526595143e-4,1.809654575597889e-2,1.817026506944104e-2,-3.155135544497688e-2,8.854765786966081e-3];
    mt19 = [-2.70566404602225e-3,9.072081989014763e-4,6.797608031085048e-4,-3.708132395431223e-4,-7.457850772355046e-5,-1.301507150247622e-3];
    mt20 = [2.69959588679965e-3,-8.840504828303794e-4,5.555858637460534e-5,-2.236279286783707e-5,6.456733767198091e-6,1.472986243782767e-5];
    mt21 = [-3.297548767579912e-5,2.582871891348531e-5,-5.354796868593621e-5,1.94342289457779e-5,-1.381619620340242e-6,3.644697108049443e-7];
    mt22 = [4.231883796523993e-6,-7.31914179169807e-7,-6.039353884043649e-6,3.671166735172252e-6,-2.791490204495906e-6,4.500318681228498e-7];
    mt23 = [3.293111524251936e-6,-1.369159179704089e-6,-1.085706780479506e-6,1.133363978308446e-6,5.897210116603036e-6,-4.530948877514295e-6];
    mt24 = [3.537233659322586e-6,-5.424788975174539e-7,5.012823405138055e-6,-2.415660532782211e-6,2.976883323723626e-6,1.824110476325179e-6];
    mt25 = [2.884210873589682e-6,-5.219786930985809e-6,4.971376335068416e-6,-9.371988579376849e-7,-3.645993820552509e-1,5.85500435115889e-1];
    mt26 = [-5.809038059472703e-1,2.86897334541382e-1,-1.115820087766582e-1,-1.460062706722942e-1,9.58791194978124e-2];
    mt27 = [-7.845772754811735e-3,5.404512456189541e-2,-1.073585336216455e-1,2.196021120334431e-1,7.032261148600782e-2];
    mt28 = [-4.089958344082099e-2,1.216020135501855e-2,-6.971374991025135e-2,3.029308173266251e-2,1.661075250481592e-2,-4.372823777791905e-4];
    mt29 = [-2.448771056945431e-2,-3.457962441712478e-3,1.167049534682351e-2,7.595243359297257e-3,-8.40727260458818e-3,1.512427423154121e-3];
    mt30 = [-1.871584050808167e-3,5.83259040533837e-4,6.890000246836112e-4,-1.859065245906191e-4,-2.618442339741163e-4,-9.402984368529219e-4];
    mt31 = [1.844215115835569e-3,-5.683104467121554e-4,5.298860713666901e-5,-2.145398783364486e-5,1.136407080476074e-5,1.293250241488135e-5];
    mt32 = [-3.518794582180465e-5,2.943669825768211e-5,-5.667834091095652e-5,1.995966659873713e-5,1.040088879687342e-6,-3.622154227479669e-7];
    mt33 = [-1.805451071278216e-6,5.455484184925187e-7,2.525286925456318e-6,-1.967603819919991e-6,1.596015270253696e-6,-2.679602602857656e-7];
    mt34 = [-1.845743483769846e-6,6.350251543472469e-7,2.52494414279991e-6,-5.845425065290343e-7,-6.356614977029584e-6,3.577093846212675e-6];
    mt35 = [-2.518579627122339e-6,3.487463570070541e-7,1.76786261591129e-6,-5.308667953443002e-7,-2.860678018742589e-6,2.071068152926695e-7];
    mt36 = [9.275231435477831e-6,-4.029719153204612e-6,2.407967023218255e-6,-2.442974490216273e-7,-5.843418698649201e-1,2.747283778763927e-1];
    mt37 = [1.414906276108596,6.777075813941426e-2,5.47838058604561e-1,1.531313734817606e-1,-2.029404896257611e-1];
    mt38 = [3.222696441025714e-2,8.932515568727145e-2,-7.426306169427307e-2,9.540973158470301e-2,1.330856181204716e-2,1.350761446954319e-2];
    mt39 = [5.772325346494919e-2,-5.857941329452624e-2,1.306092327409366e-2,-6.265686118651367e-4,4.33818651172865e-3,-1.496408459387288e-2];
    mt40 = [-4.356802989107091e-3,5.535678359583332e-3,-4.717829337002981e-5,5.898261644569856e-3,-2.395301446171417e-3,-1.037380008281677e-3];
    mt41 = [2.591677560876387e-4,6.978536565301408e-4,-3.582887903991997e-7,-4.504785257222355e-4,-5.804581534477289e-4,9.911312680272123e-4];
    mt42 = [-2.532445148368363e-4,5.226722109448133e-5,-2.110647369096428e-5,1.065396348894746e-5,1.271839444583237e-5,-3.468028852423276e-5];
    mt43 = [2.945529814259529e-5,-5.658540343260007e-5,1.989153004841845e-5,2.664161685791366e-6,-1.00186404207108e-6,-3.321900465047802e-6];
    mt44 = [1.349823162660167e-6,4.646876573970309e-6,-4.171963227573271e-6,3.415727624845222e-6,-5.614725832941004e-7,-4.0545469672308e-6];
    mt45 = [1.784650377680525e-6,2.773086618006362e-7,-1.584804676858715e-6,-4.664744105477731e-6,5.025472296352544e-6,-4.332859325767232e-6];
    mt46 = [7.447818991232361e-7,-9.063781094932684e-7,2.224133640344991e-7,1.320654753304626e-6,2.808810837674463e-7,-6.91097585494802e-6,2.072037068393604e-6];
    mt47 = [-6.788854006749521e-7,-8.635488818588362e-8,1.557774354770308e-1,-7.875940530247277e-2,-2.668308130689734e-1,5.357392671647787e-1];
    mt48 = [-1.900091375178245e-1,3.097516765775241e-1,-1.802245061378202e-1,2.741306812040095e-2,3.736262885803773e-2,-2.452378387537094e-2];
    mt49 = [3.809326671438695e-2,-3.398099049067763e-2,2.585922585115711e-2,5.282391848491035e-2,8.5006899938756e-3,-7.153076827860643e-3];
    mt50 = [-1.277829559996691e-2,7.468301624211811e-3,-8.070768682254016e-3,-4.731212322075497e-3,7.027116916772896e-4,-3.118572535730116e-3];
    mt51 = [1.012911132450378e-2,-2.492824412378216e-3,1.911278995798304e-3,-9.88125970481346e-4,-9.077884722244441e-5,1.286102869907715e-3];
    mt52 = [-1.01515868201163e-3,-7.581380264479763e-4,-1.597603977421372e-3,2.135562491168229e-3,8.74431857628556e-3,-3.840582090100452e-3];
    mt53 = [-1.210664992390249e-3,4.187100326001495e-3,-4.410198251506972e-3,-8.519593503418717e-4,-6.036672571158033e-3,6.251430189736779e-3];
    mt54 = [-2.373693648526872e-4,1.066328684299147e-4,1.717222626798273e-5,-1.161985820607036e-4,1.165079897464498e-4,2.68259543091117e-5];
    mt55 = [1.896333831279177e-4,-2.041889006388878e-4,1.158890600527161e-4,-5.12768268139633e-5,-7.755407637895099e-6,5.440236664852187e-5];
    mt56 = [-6.052619507096742e-5,-1.091028476303578e-5,-7.637357008545235e-5,7.916995249872878e-5,-1.463959258489728e-4,6.562182531897641e-5];
    mt57 = [4.275977325501953e-5,-8.375111948027461e-5,9.166017610084429e-5,3.738539986050216e-5,8.306962868792484e-5,-1.155816088235965e-4];
    mt58 = [3.036465055850314e-1,-1.603137434203871e-1,5.321572176304301e-3,2.39092389875523e-1,9.588568718240866e-1];
    mt59 = [-5.526091701382724e-2,6.09867464104082e-1,-1.569106678589582e-1,-5.739235308066453e-2,3.128139701239055e-2];
    mt60 = [-8.260317374572069e-3,-7.239972959833572e-2,5.234237382746589e-2,2.888240199532394e-2,7.159891121584665e-2,-3.192721282792944e-2];
    mt61 = [1.062821938304614e-2,-5.176431440267871e-3,9.158539158962998e-3,7.904676769109529e-3,-1.526430986808955e-2,-1.697185058760339e-2];
    mt62 = [2.202269118416011e-2,-5.505907904821996e-3,-1.047743795881001e-2,4.476939315507908e-3,5.900051991551498e-3,-5.825201783985429e-3];
    mt63 = [6.107711642862084e-4,4.28746693548744e-3,8.34987978877995e-3,-1.322380919682949e-2,-2.906744565394828e-2,1.301145125269944e-2];
    mt64 = [4.852393636832012e-3,-1.545529514817613e-2,1.677181954689432e-2,6.78144070367373e-3,1.787748321050376e-2,-2.715210692131096e-2];
    mt65 = [1.251189262515766e-3,-5.685397404297861e-4,-1.685477269936809e-4,6.763131811815047e-4,-7.229178713222715e-4,-2.667144285932205e-4];
    mt66 = [-9.447059677999236e-4,1.309337411480877e-3,-3.838649789020297e-4,1.735965251355502e-4,4.011153283105188e-5,-2.066378773554884e-4];
    mt67 = [2.773303088614513e-4,7.665111652352699e-5,2.2709483161151e-4,-3.428019605096663e-4,7.768297436095363e-4,-3.529159536706636e-4];
    mt68 = [-1.926840006471602e-4,4.50993079264632e-4,-4.393133108059863e-4,-2.641434790466516e-4,-4.419363966360603e-4,7.637399732311855e-4];
    mt69 = [-4.837270231433454e-2,4.693429630313338e-2,-1.782122262760429e-1,-6.809759940365773e-2,-3.701161748973242e-2];
    mt70 = [6.780798322309028e-1,2.387487043380946e-1,1.444707529494299e-2,-2.106715453350842e-2,-1.556255429523533e-3,1.461865909886239e-1];
    mt71 = [-8.747257150417628e-3,-1.374870943850936e-1,-7.484148727297132e-2,1.125463450230208e-1,6.6011255253836e-2];
    mt72 = [-1.086401325041653e-2,5.260790697935942e-3,-4.922189415143839e-3,-4.2358922780723e-3,1.917436249225498e-3,1.775399625928922e-2];
    mt73 = [-4.430626920577541e-2,2.667092984374549e-2,2.584897761640898e-2,-1.153629448132959e-2,-1.26714877551619e-2,1.550870036621203e-2];
    mt74 = [-3.50521065384716e-3,-1.301647833128929e-2,-2.05552449198034e-2,4.014347212308127e-2,2.784010478989635e-2,-1.30109564815765e-2];
    mt75 = [-7.10995806169632e-3,1.839627113145262e-2,-1.893349569735843e-2,-1.767778435695396e-2,-9.407784783324541e-3,4.369758889892357e-2];
    mt76 = [-2.734051414647051e-3,1.269114402944894e-3,5.761251925621583e-4,-1.700656089012951e-3,2.06024303677461e-3,1.007494143418471e-3];
    mt77 = [1.802869885519316e-3,-3.688113133149875e-3,3.758090091022605e-4,-1.77510351983398e-4,-1.098162744405741e-4,2.690498949892713e-4];
    mt78 = [-3.854765593806866e-4,-2.317098407709113e-4,-8.127914419732646e-5,5.41279573411777e-4,-1.730473747179165e-3,8.007985798059738e-4];
    mt79 = [2.992332382347504e-4,-1.011991265053235e-3,7.333369995684263e-4,8.404628172876336e-4,9.252195052890888e-4,-2.214820745687691e-3];
    mt80 = [-8.745307911975908e-2,4.018653166948424e-2,6.38864145824049e-2,-6.010633766901577e-2,-1.093160185427294e-1];
    mt81 = [1.560759090362972e-1,3.71293951371357e-1,6.28215091440208e-1,-1.433198149042145e-2,1.240459476194677e-2];
    mt82 = [-7.392535324461194e-2,-1.057735742683476e-2,1.134034567870015e-1,-3.936181166511445e-2,-9.019034298147026e-2];
    mt83 = [5.030759113361178e-2,3.208863149118665e-3,1.066239695656385e-3,-4.183665243359012e-2,2.456063173659831e-3,5.357158818609069e-2];
    mt84 = [-1.646970410804708e-2,7.663990287339727e-3,-1.985589266541402e-2,-2.972580463728368e-2,1.443183593112895e-2,9.205458657564858e-3];
    mt85 = [-2.183454806296911e-2,1.535123067976475e-2,2.529381342595623e-2,1.608545665575387e-2,-6.697004472971334e-2,1.191187582684787e-2];
    mt86 = [-4.395364976199355e-3,2.761825643080945e-3,5.712394975816661e-4,-8.541530747114069e-3,1.887149499263738e-2,-2.267018650070367e-2];
    mt87 = [-2.257355279856001e-2,3.027102779058806e-3,-1.45722924023215e-3,-1.179733591971947e-3,2.352136893535683e-3,-3.068605565070656e-3];
    mt88 = [-2.125460756945824e-3,-1.304460499345801e-3,5.943893865711988e-3,1.414795310538727e-4,-5.122912472243588e-5,1.09099158354695e-4];
    mt89 = [-3.180706843383086e-5,1.542138881021123e-5,2.655197620756998e-4,-3.482066830426838e-4,-2.678042976688708e-4,1.996472773878987e-3];
    mt90 = [-9.510515921373981e-4,-8.162570674172113e-5,1.159123529529056e-3,-1.978949905235366e-4,-1.493016257284978e-3,-9.511736108680862e-4];
    mt91 = [3.703606146217452e-3,4.596870862013336e-2,-1.124501886313701e-2,-1.204150699022877e-1,1.002062799064334e-2,2.491388189851239e-1];
    mt92 = [-9.337618138876158e-2,-3.778782309643057e-2,5.627788359704207e-1,6.955372984213398e-2,-8.319821043227076e-3,-2.260791192884091e-1];
    mt93 = [-3.301392282927394e-3,3.298235320497532e-1,-1.122359162857098e-2,-1.821688450554259e-1,7.783514361571917e-3,1.97334103814388e-2];
    mt94 = [-1.048661539177586e-2,3.044023809436101e-2,7.428783744015015e-3,-3.539142170169491e-2,-4.611787214869856e-5,1.957257588409096e-2];
    mt95 = [-1.428691290142984e-2,1.203399084017328e-2,-7.534357952463405e-3,-4.716183997189138e-4,1.53766872737242e-2,-1.668722385684858e-2];
    mt96 = [-2.842466363426004e-2,3.673149375148373e-3,6.780691075361767e-2,-4.702076226555912e-2,2.074483689193516e-2,3.55483342030854e-3];
    mt97 = [-2.343353605954639e-2,4.110699423547737e-2,3.359270244567499e-4,3.047156314928984e-2,-2.188850155859401e-2,-1.316609496754684e-3];
    mt98 = [7.116452087222886e-4,1.545805227439146e-3,-1.75768958459775e-3,2.032308741899122e-3,2.676476611843407e-3,-4.048240535429701e-4];
    mt99 = [-6.027005874792892e-3,-6.014076690347304e-4,2.696462206113578e-4,4.169864441144553e-6,-3.086601957619658e-4,6.182968339575037e-4];
    mt100 = [-3.307374881598654e-6,3.702741503408926e-4,-2.704187947104768e-4,-1.081494520317039e-3,5.531936269286985e-4,-2.702182926124177e-4];
    mt101 = [-6.099906517285368e-4,-1.09087085646441e-3,1.650258784523307e-3,3.557861749268566e-4,-3.898434113440164e-3,-2.401806647083848e-2];
    mt102 = [-1.528456709086159e-2,2.277082724768517e-1,4.183954526659336e-2,-4.170876622594765e-1,-2.130242008915423e-2];
    mt103 = [2.007687073421834e-1,-1.631415011219778e-1,-3.897906374081648e-2,1.405362067730507e-2,4.039272700286246e-2];
    mt104 = [-1.204419894649187e-2,-1.035853371534702e-1,5.251933077525049e-2,4.206369512613188e-3,-8.553784681100628e-2];
    mt105 = [-1.883971697025454e-2,9.398020549772634e-3,-8.283396368658571e-3,-1.027713642226229e-2,6.630149252082612e-3,2.08329487636472e-2];
    mt106 = [-3.059385664856707e-2,2.615464700488442e-2,8.432599890695919e-3,-2.818241411993311e-3,6.546745413776895e-3,-1.821826751388475e-3];
    mt107 = [-1.005873236635501e-2,1.852771875911753e-2,-8.92116394931663e-3,-4.51820777928667e-2,2.834342908811432e-2];
    mt108 = [-1.431418740721278e-2,-2.182895621756264e-3,2.253342011145291e-2,-3.38303912353947e-2,-2.853994927744252e-2];
    mt109 = [1.334869517903585e-2,3.588478672228797e-2,-6.266902637519811e-4,2.285145634718811e-4,-1.171417366686861e-3,2.895458965077058e-4];
    mt110 = [8.917764267118046e-4,-1.718938698603533e-3,6.962205864918686e-4,4.171134302642487e-3,3.630267421931991e-4,-1.873576330598501e-4];
    mt111 = [-1.50047864829854e-4,3.656373275913731e-4,-7.160318134192292e-4,-4.248870016190383e-4,3.26862372044254e-4,4.046503763775639e-4];
    mt112 = [7.536190415922844e-5,-6.069766331592087e-5,6.942172126017606e-6,4.245876948032465e-5,1.730283501472471e-3,-1.152926944370197e-3];
    mt113 = [1.246057142923107e-4,2.682475664243514e-3,-4.050058807515699e-4,7.975927014411752e-3,-4.072547745157325e-2,-2.650352873408943e-2];
    mt114 = [1.363524508858941e-1,3.338265524402222e-2,-9.595130409112441e-2,5.950203064743431e-2,1.204358880900257e-2,-1.751517164502629e-2];
    mt115 = [1.271793507789728e-1,1.974334214728427e-2,-1.436094885338519e-1,-2.754510756520812e-2,7.586458552787288e-2];
    mt116 = [1.070359272529913e-2,-1.567625063657356e-2,1.816126909629524e-3,4.264018129185302e-2,3.611367588188309e-3,-7.117310039864189e-2];
    mt117 = [-4.133461018316991e-3,3.76848719462476e-2,-2.100018905504086e-2,-7.248517572207595e-3,4.315475975993333e-3,-1.434782105743885e-2];
    mt118 = [-4.283120421038701e-3,2.63479781033512e-2,-2.121101556934588e-3,-7.524065049385301e-3,2.361022772447176e-2,2.180066196335822e-2];
    mt119 = [-6.534732484116744e-3,-7.392733620672055e-3,-7.128424809399764e-4,-3.501271563926551e-3,3.429300104243161e-2,-4.882806578918193e-2];
    mt120 = [-1.50596717893706e-2,5.607604663945645e-4,-3.468300113818716e-4,2.365894457849571e-4,7.661674145811988e-4,-3.20837682882798e-3];
    mt121 = [-3.89728876591711e-4,1.727417021905674e-3,-2.522937675088609e-3,2.57351818923206e-4,-7.599131788304326e-5,1.141530589301552e-4];
    mt122 = [-9.75695646653349e-5,1.709003510074351e-4,5.243310821225605e-4,-7.818176103553824e-4,-1.430556964196357e-4,-2.288486103429953e-4];
    mt123 = [3.01579222714846e-5,1.161589113330643e-3,-1.664102832770666e-4,-8.526879273108953e-4,4.306404620429402e-4,1.206028739578741e-4];
    mt124 = [-1.31133756361025e-3,-8.221097311847905e-3,2.24948264595914e-3,-7.782781425462211e-3,7.105412313102793e-3,-3.600276332681467e-2];
    mt125 = [-9.867585475358077e-3,2.798540954380206e-2,-1.601979169432427e-2,-4.080403355571926e-2,3.163828445367605e-3,1.345010835514535e-1];
    mt126 = [9.260892975690341e-3,-2.139098359856896e-1,-8.308792177786871e-4,9.058152964218211e-2,-4.702208577250475e-2];
    mt127 = [2.473896468052687e-2,-1.202567400792914e-2,1.186166022504137e-2,1.124609626268011e-2,4.749158371637619e-4,-2.380987313346207e-2];
    mt128 = [2.488539912799442e-2,-3.508302345072522e-3,-2.054083990478769e-3,-5.825395188448153e-4,1.13207773397171e-2,2.379642323319201e-3];
    mt129 = [-1.765104395871099e-2,-5.406763988921559e-3,1.538086204387184e-2,-1.07215575716398e-2,-7.515874572087376e-4,4.794944473180153e-4];
    mt130 = [1.127638176670251e-2,-5.560985850742248e-3,9.913892830965287e-3,9.675869024396473e-3,-1.747547775829393e-2,1.667213366304006e-2];
    mt131 = [1.667668195043459e-4,9.03688183279397e-5,8.214230018362737e-5,-8.411560524790137e-4,3.294625458275359e-3,1.31673939573025e-3];
    mt132 = [-2.960851202291939e-3,1.46540188255511e-3,4.453145052981102e-5,-1.554819118745058e-5,2.135003714468711e-4,-8.580049249849493e-5];
    mt133 = [1.993670807876718e-4,1.63657060214563e-4,-3.193381010242055e-4,2.443941419805057e-4,5.090296757277146e-4,-9.864943312012773e-5];
    mt134 = [-2.158621559818065e-3,4.508741043134848e-4,1.720041507715224e-4,-3.241843501425658e-4,-3.559924452310608e-5,4.136859998587163e-4];
    mt135 = [7.706628741162949e-3,1.226393585043829e-3,-4.065344108473701e-2,-6.739930963430219e-3,8.765510761226422e-2,-5.429367337169651e-3];
    mt136 = [-2.271417924231228e-2,1.081895091792072e-2,3.202799552583988e-2,-1.256675188593808e-2,-1.165174457301051e-2,1.440849526123545e-2];
    mt137 = [-5.023927304286893e-2,2.034589903809585e-2,-3.198124325200028e-2,1.21003586938578e-2,6.897269887634756e-2];
    mt138 = [-3.081056302590946e-2,5.394575627571113e-2,1.478452080934135e-2,-4.103354284880743e-2,2.378704619634806e-2,-4.99553401989411e-2];
    mt139 = [1.84903448762981e-2,6.665099247081835e-4,-5.204428775988141e-4,9.199924696108936e-3,-1.945069511084353e-3,-1.681097060371926e-3];
    mt140 = [5.603454642055482e-3,-7.26841504223506e-3,2.580498071109898e-3,-3.708566163810063e-2,1.308281549954089e-2,1.171453420980623e-2];
    mt141 = [-6.606009630803499e-3,2.198033500662896e-3,-2.452981712610935e-2,5.331278088220026e-2,-2.050328374812438e-2,8.702479306392977e-4];
    mt142 = [-4.200941405425471e-4,1.916370127548403e-4,5.471368172621063e-4,-1.433667123518629e-3,-3.098949501696728e-4,4.709079267691432e-4];
    mt143 = [-1.547443768630711e-4,-4.325902316662777e-4,1.501267315846127e-4,-1.216530103426741e-4,4.979159210065333e-5,-3.903842001995936e-4];
    mt144 = [-3.891410025251949e-4,8.615882042470557e-4,-3.158533419768342e-4,4.372001073797985e-4,-2.957829948465651e-4,1.864558707227052e-3];
    mt145 = [7.393807562737986e-5,-1.984284564472035e-3,7.200075182015132e-4,-5.534663496788862e-4,8.106628231552405e-5];
    BCU_g2f = reshape([mt15,mt16,mt17,mt18,mt19,mt20,mt21,mt22,mt23,mt24,mt25,mt26,mt27,mt28,mt29,mt30,mt31,mt32,mt33,mt34,mt35,mt36,mt37,mt38,mt39,mt40,mt41,mt42,mt43,mt44,mt45,mt46,mt47,mt48,mt49,mt50,mt51,mt52,mt53,mt54,mt55,mt56,mt57,mt58,mt59,mt60,mt61,mt62,mt63,mt64,mt65,mt66,mt67,mt68,mt69,mt70,mt71,mt72,mt73,mt74,mt75,mt76,mt77,mt78,mt79,mt80,mt81,mt82,mt83,mt84,mt85,mt86,mt87,mt88,mt89,mt90,mt91,mt92,mt93,mt94,mt95,mt96,mt97,mt98,mt99,mt100,mt101,mt102,mt103,mt104,mt105,mt106,mt107,mt108,mt109,mt110,mt111,mt112,mt113,mt114,mt115,mt116,mt117,mt118,mt119,mt120,mt121,mt122,mt123,mt124,mt125,mt126,mt127,mt128,mt129,mt130,mt131,mt132,mt133,mt134,mt135,mt136,mt137,mt138,mt139,mt140,mt141,mt142,mt143,mt144,mt145],8,96);
end
if nargout > 2
    mt146 = [2.948906761778786e-1,1.525720623897707,2.57452876984127e-1,1.798113701499118,4.127080577601411e-1,1.278484623015873];
    mt147 = [9.232955798059965e-1,1.009333860859158];
    HU = reshape([mt146,mt147],8,1);
end
if nargout > 3
    M = [2.0;2.0./3.0;2.0./5.0;2.0./7.0;2.0./9.0;2.0./1.1e+1;2.0./1.3e+1;2.0./1.5e+1];
end