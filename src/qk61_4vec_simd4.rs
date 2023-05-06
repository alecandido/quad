use std::iter::zip;
use std::simd::{f64x4,SimdFloat};
use simba::simd::Simd;
use crate::funct_vector::FnVec4;
use crate::qk::*;

pub struct Qk61Vec4Simd4 {}
pub const XGK: [f64;61] = [-0.999484410050490637571325895705811, -0.996893484074649540271630050918695,
    -0.991630996870404594858628366109486, -0.983668123279747209970032581605663,
    -0.973116322501126268374693868423707, -0.960021864968307512216871025581798,
    -0.944374444748559979415831324037439, -0.926200047429274325879324277080474,
    -0.905573307699907798546522558925958, -0.882560535792052681543116462530226,
    -0.857205233546061098958658510658944, -0.829565762382768397442898119732502,
    -0.799727835821839083013668942322683, -0.767777432104826194917977340974503,
    -0.733790062453226804726171131369528, -0.697850494793315796932292388026640,
    -0.660061064126626961370053668149271, -0.620526182989242861140477556431189,
    -0.579345235826361691756024932172540, -0.536624148142019899264169793311073,
    -0.492480467861778574993693061207709, -0.447033769538089176780609900322854,
    -0.400401254830394392535476211542661, -0.352704725530878113471037207089374,
    -0.304073202273625077372677107199257, -0.254636926167889846439805129817805,
    -0.204525116682309891438957671002025, -0.153869913608583546963794672743256,
    -0.102806937966737030147096751318001, -0.051471842555317695833025213166723,
    0.000000000000000000000000000000000, 0.051471842555317695833025213166723,
    0.102806937966737030147096751318001, 0.153869913608583546963794672743256,
    0.204525116682309891438957671002025, 0.254636926167889846439805129817805,
    0.304073202273625077372677107199257, 0.352704725530878113471037207089374,
    0.400401254830394392535476211542661, 0.447033769538089176780609900322854,
    0.492480467861778574993693061207709, 0.536624148142019899264169793311073,
    0.579345235826361691756024932172540, 0.620526182989242861140477556431189,
    0.660061064126626961370053668149271, 0.697850494793315796932292388026640,
    0.733790062453226804726171131369528, 0.767777432104826194917977340974503,
    0.799727835821839083013668942322683, 0.829565762382768397442898119732502,
    0.857205233546061098958658510658944, 0.882560535792052681543116462530226,
    0.905573307699907798546522558925958, 0.926200047429274325879324277080474,
    0.944374444748559979415831324037439, 0.960021864968307512216871025581798,
    0.973116322501126268374693868423707, 0.983668123279747209970032581605663,
    0.991630996870404594858628366109486, 0.996893484074649540271630050918695,
    0.999484410050490637571325895705811];


pub const WGK1: [f64x4;31] = [f64x4::from_array([0.001389013698677007624551591226760,0.001389013698677007624551591226760,0.001389013698677007624551591226760,0.001389013698677007624551591226760]),
    f64x4::from_array([0.006630703915931292173319826369750,0.006630703915931292173319826369750,0.006630703915931292173319826369750,0.006630703915931292173319826369750]),
    f64x4::from_array([0.011823015253496341742232898853251,0.011823015253496341742232898853251,0.011823015253496341742232898853251,0.011823015253496341742232898853251]),
    f64x4::from_array([0.016920889189053272627572289420322,0.016920889189053272627572289420322,0.016920889189053272627572289420322,0.016920889189053272627572289420322]),
    f64x4::from_array([0.021828035821609192297167485738339,0.021828035821609192297167485738339,0.021828035821609192297167485738339,0.021828035821609192297167485738339]),
    f64x4::from_array([0.026509954882333101610601709335075,0.026509954882333101610601709335075,0.026509954882333101610601709335075,0.026509954882333101610601709335075]),
    f64x4::from_array([0.030907257562387762472884252943092,0.030907257562387762472884252943092,0.030907257562387762472884252943092,0.030907257562387762472884252943092]),
    f64x4::from_array([0.034979338028060024137499670731468,0.034979338028060024137499670731468,0.034979338028060024137499670731468,0.034979338028060024137499670731468]),
    f64x4::from_array([0.038678945624727592950348651532281,0.038678945624727592950348651532281,0.038678945624727592950348651532281,0.038678945624727592950348651532281]),
    f64x4::from_array([0.041969810215164246147147541285970,0.041969810215164246147147541285970,0.041969810215164246147147541285970,0.041969810215164246147147541285970]),
    f64x4::from_array([0.044814800133162663192355551616723,0.044814800133162663192355551616723,0.044814800133162663192355551616723,0.044814800133162663192355551616723]),
    f64x4::from_array([0.047185546569299153945261478181099,0.047185546569299153945261478181099,0.047185546569299153945261478181099,0.047185546569299153945261478181099]),
    f64x4::from_array([0.049055434555029778887528165367238,0.049055434555029778887528165367238,0.049055434555029778887528165367238,0.049055434555029778887528165367238]),
    f64x4::from_array([0.050405921402782346840893085653585,0.050405921402782346840893085653585,0.050405921402782346840893085653585,0.050405921402782346840893085653585]),
    f64x4::from_array([0.051221547849258772170656282604944,0.051221547849258772170656282604944,0.051221547849258772170656282604944,0.051221547849258772170656282604944]),
    f64x4::from_array([0.051494729429451567558340433647099,0.051494729429451567558340433647099,0.051494729429451567558340433647099,0.051494729429451567558340433647099]),
    f64x4::from_array([0.051221547849258772170656282604944,0.051221547849258772170656282604944,0.051221547849258772170656282604944,0.051221547849258772170656282604944]),
    f64x4::from_array([0.050405921402782346840893085653585,0.050405921402782346840893085653585,0.050405921402782346840893085653585,0.050405921402782346840893085653585]),
    f64x4::from_array([0.049055434555029778887528165367238,0.049055434555029778887528165367238,0.049055434555029778887528165367238,0.049055434555029778887528165367238]),
    f64x4::from_array([0.047185546569299153945261478181099,0.047185546569299153945261478181099,0.047185546569299153945261478181099,0.047185546569299153945261478181099]),
    f64x4::from_array([0.044814800133162663192355551616723,0.044814800133162663192355551616723,0.044814800133162663192355551616723,0.044814800133162663192355551616723]),
    f64x4::from_array([0.041969810215164246147147541285970,0.041969810215164246147147541285970,0.041969810215164246147147541285970,0.041969810215164246147147541285970]),
    f64x4::from_array([0.038678945624727592950348651532281,0.038678945624727592950348651532281,0.038678945624727592950348651532281,0.038678945624727592950348651532281]),
    f64x4::from_array([0.034979338028060024137499670731468,0.034979338028060024137499670731468,0.034979338028060024137499670731468,0.034979338028060024137499670731468]),
    f64x4::from_array([0.030907257562387762472884252943092,0.030907257562387762472884252943092,0.030907257562387762472884252943092,0.030907257562387762472884252943092]),
    f64x4::from_array([0.026509954882333101610601709335075,0.026509954882333101610601709335075,0.026509954882333101610601709335075,0.026509954882333101610601709335075]),
    f64x4::from_array([0.021828035821609192297167485738339,0.021828035821609192297167485738339,0.021828035821609192297167485738339,0.021828035821609192297167485738339]),
    f64x4::from_array([0.016920889189053272627572289420322,0.016920889189053272627572289420322,0.016920889189053272627572289420322,0.016920889189053272627572289420322]),
    f64x4::from_array([0.011823015253496341742232898853251,0.011823015253496341742232898853251,0.011823015253496341742232898853251,0.011823015253496341742232898853251]),
    f64x4::from_array([0.006630703915931292173319826369750,0.006630703915931292173319826369750,0.006630703915931292173319826369750,0.006630703915931292173319826369750]),
    f64x4::from_array([0.001389013698677007624551591226760,0.001389013698677007624551591226760,0.001389013698677007624551591226760,0.001389013698677007624551591226760])];

pub const WGK2: [f64x4;31] = [f64x4::from_array([0.003890461127099884051267201844516,0.003890461127099884051267201844516,0.003890461127099884051267201844516,0.003890461127099884051267201844516]),
    f64x4::from_array([0.009273279659517763428441146892024,0.009273279659517763428441146892024,0.009273279659517763428441146892024,0.009273279659517763428441146892024]),
    f64x4::from_array([0.014369729507045804812451432443580,0.014369729507045804812451432443580,0.014369729507045804812451432443580,0.014369729507045804812451432443580]),
    f64x4::from_array([0.019414141193942381173408951050128,0.019414141193942381173408951050128,0.019414141193942381173408951050128,0.019414141193942381173408951050128]),
    f64x4::from_array([0.024191162078080601365686370725232,0.024191162078080601365686370725232,0.024191162078080601365686370725232,0.024191162078080601365686370725232]),
    f64x4::from_array([0.028754048765041292843978785354334,0.028754048765041292843978785354334,0.028754048765041292843978785354334,0.028754048765041292843978785354334]),
    f64x4::from_array([0.032981447057483726031814191016854,0.032981447057483726031814191016854,0.032981447057483726031814191016854,0.032981447057483726031814191016854]),
    f64x4::from_array([0.036882364651821229223911065617136,0.036882364651821229223911065617136,0.036882364651821229223911065617136,0.036882364651821229223911065617136]),
    f64x4::from_array([0.040374538951535959111995279752468,0.040374538951535959111995279752468,0.040374538951535959111995279752468,0.040374538951535959111995279752468]),
    f64x4::from_array([0.043452539701356069316831728117073,0.043452539701356069316831728117073,0.043452539701356069316831728117073,0.043452539701356069316831728117073]),
    f64x4::from_array([0.046059238271006988116271735559374,0.046059238271006988116271735559374,0.046059238271006988116271735559374,0.046059238271006988116271735559374]),
    f64x4::from_array([0.048185861757087129140779492298305,0.048185861757087129140779492298305,0.048185861757087129140779492298305,0.048185861757087129140779492298305]),
    f64x4::from_array([0.049795683427074206357811569379942,0.049795683427074206357811569379942,0.049795683427074206357811569379942,0.049795683427074206357811569379942]),
    f64x4::from_array([0.050881795898749606492297473049805,0.050881795898749606492297473049805,0.050881795898749606492297473049805,0.050881795898749606492297473049805]),
    f64x4::from_array([0.051426128537459025933862879215781,0.051426128537459025933862879215781,0.051426128537459025933862879215781,0.051426128537459025933862879215781]),
    f64x4::from_array([0.051426128537459025933862879215781,0.051426128537459025933862879215781,0.051426128537459025933862879215781,0.051426128537459025933862879215781]),
    f64x4::from_array([0.050881795898749606492297473049805,0.050881795898749606492297473049805,0.050881795898749606492297473049805,0.050881795898749606492297473049805]),
    f64x4::from_array([0.049795683427074206357811569379942,0.049795683427074206357811569379942,0.049795683427074206357811569379942,0.049795683427074206357811569379942]),
    f64x4::from_array([0.048185861757087129140779492298305,0.048185861757087129140779492298305,0.048185861757087129140779492298305,0.048185861757087129140779492298305]),
    f64x4::from_array([0.046059238271006988116271735559374,0.046059238271006988116271735559374,0.046059238271006988116271735559374,0.046059238271006988116271735559374]),
    f64x4::from_array([0.043452539701356069316831728117073,0.043452539701356069316831728117073,0.043452539701356069316831728117073,0.043452539701356069316831728117073]),
    f64x4::from_array([0.040374538951535959111995279752468,0.040374538951535959111995279752468,0.040374538951535959111995279752468,0.040374538951535959111995279752468]),
    f64x4::from_array([0.036882364651821229223911065617136,0.036882364651821229223911065617136,0.036882364651821229223911065617136,0.036882364651821229223911065617136]),
    f64x4::from_array([0.032981447057483726031814191016854,0.032981447057483726031814191016854,0.032981447057483726031814191016854,0.032981447057483726031814191016854]),
    f64x4::from_array([0.028754048765041292843978785354334,0.028754048765041292843978785354334,0.028754048765041292843978785354334,0.028754048765041292843978785354334]),
    f64x4::from_array([0.024191162078080601365686370725232,0.024191162078080601365686370725232,0.024191162078080601365686370725232,0.024191162078080601365686370725232]),
    f64x4::from_array([0.019414141193942381173408951050128,0.019414141193942381173408951050128,0.019414141193942381173408951050128,0.019414141193942381173408951050128]),
    f64x4::from_array([0.014369729507045804812451432443580,0.014369729507045804812451432443580,0.014369729507045804812451432443580,0.014369729507045804812451432443580]),
    f64x4::from_array([0.009273279659517763428441146892024,0.009273279659517763428441146892024,0.009273279659517763428441146892024,0.009273279659517763428441146892024]),
    f64x4::from_array([0.003890461127099884051267201844516,0.003890461127099884051267201844516,0.003890461127099884051267201844516,0.003890461127099884051267201844516]),
    f64x4::from_array([0.000000000000000000000000000000000,0.000000000000000000000000000000000,0.000000000000000000000000000000000,0.000000000000000000000000000000000])];

pub const WG: [f64x4;31] = [f64x4::from_array([0.007968192496166605615465883474674,0.007968192496166605615465883474674,0.007968192496166605615465883474674,0.007968192496166605615465883474674]),
    f64x4::from_array([0.018466468311090959142302131912047,0.018466468311090959142302131912047,0.018466468311090959142302131912047,0.018466468311090959142302131912047]),
    f64x4::from_array([0.028784707883323369349719179611292,0.028784707883323369349719179611292,0.028784707883323369349719179611292,0.028784707883323369349719179611292]),
    f64x4::from_array([0.038799192569627049596801936446348,0.038799192569627049596801936446348,0.038799192569627049596801936446348,0.038799192569627049596801936446348]),
    f64x4::from_array([0.048402672830594052902938140422808,0.048402672830594052902938140422808,0.048402672830594052902938140422808,0.048402672830594052902938140422808]),
    f64x4::from_array([0.057493156217619066481721689402056,0.057493156217619066481721689402056,0.057493156217619066481721689402056,0.057493156217619066481721689402056]),
    f64x4::from_array([0.065974229882180495128128515115962,0.065974229882180495128128515115962,0.065974229882180495128128515115962,0.065974229882180495128128515115962]),
    f64x4::from_array([0.073755974737705206268243850022191,0.073755974737705206268243850022191,0.073755974737705206268243850022191,0.073755974737705206268243850022191]),
    f64x4::from_array([0.080755895229420215354694938460530,0.080755895229420215354694938460530,0.080755895229420215354694938460530,0.080755895229420215354694938460530]),
    f64x4::from_array([0.086899787201082979802387530715126,0.086899787201082979802387530715126,0.086899787201082979802387530715126,0.086899787201082979802387530715126]),
    f64x4::from_array([0.092122522237786128717632707087619,0.092122522237786128717632707087619,0.092122522237786128717632707087619,0.092122522237786128717632707087619]),
    f64x4::from_array([0.096368737174644259639468626351810,0.096368737174644259639468626351810,0.096368737174644259639468626351810,0.096368737174644259639468626351810]),
    f64x4::from_array([0.099593420586795267062780282103569,0.099593420586795267062780282103569,0.099593420586795267062780282103569,0.099593420586795267062780282103569]),
    f64x4::from_array([0.101762389748405504596428952168554,0.101762389748405504596428952168554,0.101762389748405504596428952168554,0.101762389748405504596428952168554]),
    f64x4::from_array([0.102852652893558840341285636705415,0.102852652893558840341285636705415,0.102852652893558840341285636705415,0.102852652893558840341285636705415]),
    f64x4::from_array([0.102852652893558840341285636705415,0.102852652893558840341285636705415,0.102852652893558840341285636705415,0.102852652893558840341285636705415]),
    f64x4::from_array([0.101762389748405504596428952168554,0.101762389748405504596428952168554,0.101762389748405504596428952168554,0.101762389748405504596428952168554]),
    f64x4::from_array([0.099593420586795267062780282103569,0.099593420586795267062780282103569,0.099593420586795267062780282103569,0.099593420586795267062780282103569]),
    f64x4::from_array([0.096368737174644259639468626351810,0.096368737174644259639468626351810,0.096368737174644259639468626351810,0.096368737174644259639468626351810]),
    f64x4::from_array([0.092122522237786128717632707087619,0.092122522237786128717632707087619,0.092122522237786128717632707087619,0.092122522237786128717632707087619]),
    f64x4::from_array([0.086899787201082979802387530715126,0.086899787201082979802387530715126,0.086899787201082979802387530715126,0.086899787201082979802387530715126]),
    f64x4::from_array([0.080755895229420215354694938460530,0.080755895229420215354694938460530,0.080755895229420215354694938460530,0.080755895229420215354694938460530]),
    f64x4::from_array([0.073755974737705206268243850022191,0.073755974737705206268243850022191,0.073755974737705206268243850022191,0.073755974737705206268243850022191]),
    f64x4::from_array([0.065974229882180495128128515115962,0.065974229882180495128128515115962,0.065974229882180495128128515115962,0.065974229882180495128128515115962]),
    f64x4::from_array([0.057493156217619066481721689402056,0.057493156217619066481721689402056,0.057493156217619066481721689402056,0.057493156217619066481721689402056]),
    f64x4::from_array([0.048402672830594052902938140422808,0.048402672830594052902938140422808,0.048402672830594052902938140422808,0.048402672830594052902938140422808]),
    f64x4::from_array([0.038799192569627049596801936446348,0.038799192569627049596801936446348,0.038799192569627049596801936446348,0.038799192569627049596801936446348]),
    f64x4::from_array([0.028784707883323369349719179611292,0.028784707883323369349719179611292,0.028784707883323369349719179611292,0.028784707883323369349719179611292]),
    f64x4::from_array([0.018466468311090959142302131912047,0.018466468311090959142302131912047,0.018466468311090959142302131912047,0.018466468311090959142302131912047]),
    f64x4::from_array([0.007968192496166605615465883474674,0.007968192496166605615465883474674,0.007968192496166605615465883474674,0.007968192496166605615465883474674]),
    f64x4::from_array([0.0,0.0,0.0,0.0])];



impl Qk61Vec4Simd4 {
    pub(crate) fn integrate(&self, fun: &FnVec4, a: f64, b: f64, ) -> (f64x4,f64x4,f64x4,f64x4 ) {
        let hlgth = f64x4::splat(0.5 * (b - a));
        let dhlgth = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);



        let mut resg =  f64x4::splat(0.0);
        let mut resasc =  f64x4::splat(0.0);


        let f = fvec_simd2(fun,centr,hlgth[0]);

        /*
        let resk = f.0[0] * WGK1[0] + f.1[0] * WGK2[0] + f.0[1] * WGK1[1] + f.1[1] * WGK2[1] +
            f.0[2] * WGK1[2] + f.1[2] * WGK2[2] + f.0[3] * WGK1[3] + f.1[3] * WGK2[3] +
            f.0[4] * WGK1[4] + f.1[4] * WGK2[4] + f.0[5] * WGK1[5] + f.1[5] * WGK2[5] +
            f.0[6] * WGK1[6] + f.1[6] * WGK2[6] + f.0[7] * WGK1[7] + f.1[7] * WGK2[7] +
            f.0[8] * WGK1[8] + f.1[8] * WGK2[8] + f.0[9] * WGK1[9] + f.1[9] * WGK2[9] +
            f.0[10] * WGK1[10] + f.1[10] * WGK2[10] + f.0[11] * WGK1[11] + f.1[11] * WGK2[11] +
            f.0[12] * WGK1[12] + f.1[12] * WGK2[12] + f.0[13] * WGK1[13] + f.1[13] * WGK2[13] +
            f.0[14] * WGK1[14] + f.1[14] * WGK2[14] + f.0[15] * WGK1[15] + f.1[15] * WGK2[15] +
            f.0[16] * WGK1[16] + f.1[16] * WGK2[16] + f.0[17] * WGK1[17] + f.1[17] * WGK2[17] +
            f.0[18] * WGK1[18] + f.1[18] * WGK2[18] + f.0[19] * WGK1[19] + f.1[19] * WGK2[19] +
            f.0[20] * WGK1[20] + f.1[20] * WGK2[20] + f.0[21] * WGK1[21] + f.1[21] * WGK2[21] +
            f.0[22] * WGK1[22] + f.1[22] * WGK2[22] + f.0[23] * WGK1[23] + f.1[23] * WGK2[23] +
            f.0[24] * WGK1[24] + f.1[24] * WGK2[24] + f.0[25] * WGK1[25] + f.1[25] * WGK2[25] +
            f.0[26] * WGK1[26] + f.1[26] * WGK2[26] + f.0[27] * WGK1[27] + f.1[27] * WGK2[27] +
            f.0[28] * WGK1[28] + f.1[28] * WGK2[28] + f.0[29] * WGK1[29] + f.1[29] * WGK2[29] +
            f.0[30] * WGK1[30] + f.1[30] * WGK2[30];
         */


        let resk = ( f.0[0] + f.0[30] )* WGK1[0] + ( f.1[0] + f.1[29]) * WGK2[0] +
            ( f.0[1] + f.0[29] ) * WGK1[1] + ( f.1[1] + f.1[28] ) * WGK2[1] +
            ( f.0[2] + f.0[28] ) * WGK1[2] + ( f.1[2] + f.1[27] ) * WGK2[2] +
            ( f.0[3] + f.0[27] ) * WGK1[3] + ( f.1[3] + f.1[26] ) * WGK2[3] +
            ( f.0[4] + f.0[26] ) * WGK1[4] + ( f.1[4] + f.1[25] ) * WGK2[4] +
            ( f.0[5] + f.0[25] ) * WGK1[5] + ( f.1[5] + f.1[24] ) * WGK2[5] +
            ( f.0[6] + f.0[24] ) * WGK1[6] + ( f.1[6] + f.1[23] ) * WGK2[6] +
            ( f.0[7] + f.0[23] ) * WGK1[7] + ( f.1[7] + f.1[22] ) * WGK2[7] +
            ( f.0[8] + f.0[22] ) * WGK1[8] + ( f.1[8] + f.1[21] ) * WGK2[8] +
            ( f.0[9] + f.0[21] ) * WGK1[9] + ( f.1[9] + f.1[20] ) * WGK2[9] +
            ( f.0[10] + f.0[20] ) * WGK1[10] + ( f.1[10] + f.1[19] ) * WGK2[10] +
            ( f.0[11] + f.0[19] ) * WGK1[11] + ( f.1[11] + f.1[18] )  * WGK2[11] +
            ( f.0[12] + f.0[18] ) * WGK1[12] + ( f.1[12] + f.1[17] )  * WGK2[12] +
            ( f.0[13] + f.0[17] ) * WGK1[13] + ( f.1[13] + f.1[16] )  * WGK2[13] +
            ( f.0[14] + f.0[16] ) * WGK1[14] + ( f.1[14] + f.1[15] )  * WGK2[14] +
            f.0[15] * WGK1[15] ;

        let mut resabs = ( f.0[0].abs() + f.0[30].abs() )* WGK1[0] + ( f.1[0].abs() + f.1[29].abs()) * WGK2[0] +
            ( f.0[1].abs() + f.0[29].abs() ) * WGK1[1] + ( f.1[1].abs() + f.1[28].abs() ) * WGK2[1] +
            ( f.0[2].abs() + f.0[28].abs() ) * WGK1[2] + ( f.1[2].abs() + f.1[27].abs() ) * WGK2[2] +
            ( f.0[3].abs() + f.0[27].abs() ) * WGK1[3] + ( f.1[3].abs() + f.1[26].abs() ) * WGK2[3] +
            ( f.0[4].abs() + f.0[26].abs() ) * WGK1[4] + ( f.1[4].abs() + f.1[25].abs() ) * WGK2[4] +
            ( f.0[5].abs() + f.0[25].abs() ) * WGK1[5] + ( f.1[5].abs() + f.1[24].abs() ) * WGK2[5] +
            ( f.0[6].abs() + f.0[24].abs() ) * WGK1[6] + ( f.1[6].abs() + f.1[23].abs() ) * WGK2[6] +
            ( f.0[7].abs() + f.0[23].abs() ) * WGK1[7] + ( f.1[7].abs() + f.1[22].abs() ) * WGK2[7] +
            ( f.0[8].abs() + f.0[22].abs() ) * WGK1[8] + ( f.1[8].abs() + f.1[21].abs() ) * WGK2[8] +
            ( f.0[9].abs() + f.0[21].abs() ) * WGK1[9] + ( f.1[9].abs() + f.1[20].abs() ) * WGK2[9] +
            ( f.0[10].abs() + f.0[20].abs() ) * WGK1[10] + ( f.1[10].abs() + f.1[19].abs() ) * WGK2[10] +
            ( f.0[11].abs() + f.0[19].abs() ) * WGK1[11] + ( f.1[11].abs() + f.1[18].abs() )  * WGK2[11] +
            ( f.0[12].abs() + f.0[18].abs() ) * WGK1[12] + ( f.1[12].abs() + f.1[17].abs() )  * WGK2[12] +
            ( f.0[13].abs() + f.0[17].abs() ) * WGK1[13] + ( f.1[13].abs() + f.1[16].abs() )  * WGK2[13] +
            ( f.0[14].abs() + f.0[16].abs() ) * WGK1[14] + ( f.1[14].abs() + f.1[15].abs() )  * WGK2[14] +
            f.0[15].abs() * WGK1[15] ;



        /*
                for k in 0..4 {
                    let (fv1,fv2) = fvec_simd(&fun.components[k], centr, hlgth[0]);
                    resk[k] = (fv1 * WGK1).reduce_sum() + (fv2 * WGK2).reduce_sum();
                    let reskh = resk[k] * 0.5;
                    let reskhs = Simd::from_array([reskh;32]);
                    resabs[k] = (fv1.abs() * WGK1).reduce_sum() + (fv2.abs() * WGK2).reduce_sum();
                    resg[k] = ( fv2 * WG).reduce_sum();
                    resasc[k] = ( WGK1 * ( fv1 - reskhs).abs()).reduce_sum() +
                        ( WGK2 * ( fv2 - reskhs).abs()).reduce_sum();
                }


         */
        let result = resk * hlgth;
        resabs = resabs * dhlgth;
        resasc = resasc * dhlgth;
        let mut abserr = ((resk - resg) * hlgth).abs();
        for k in 0..4 {
            if (resasc[k], abserr[k]) != (0.0, 0.0) {
                abserr[k] = resasc[k] * 1.0_f64.min((200.0 * abserr[k] / resasc[k]).powf(1.5));
            }

            if resabs[k] > UFLOW / (50.0 * EPMACH) {
                abserr[k] = abserr[k].max((EPMACH * 50.0) * resabs[k]);
            }
        }

        (result, abserr, resabs, resasc)
    }
}
pub fn fvec_simd2(fun: &FnVec4, centr : f64, hlgth : f64) -> ([f64x4;31],[f64x4;31]){
    ([f64x4::from_array([fun.components[0](centr + hlgth * XGK[0]),
                 fun.components[1](centr + hlgth * XGK[0]), fun.components[2](centr + hlgth * XGK[0]),
                 fun.components[3](centr + hlgth * XGK[0])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[2]),
                    fun.components[1](centr + hlgth * XGK[2]), fun.components[2](centr + hlgth * XGK[2]),
                    fun.components[3](centr + hlgth * XGK[2])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[4]),
                    fun.components[1](centr + hlgth * XGK[4]), fun.components[2](centr + hlgth * XGK[4]),
                    fun.components[3](centr + hlgth * XGK[4])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[6]),
                    fun.components[1](centr + hlgth * XGK[6]), fun.components[2](centr + hlgth * XGK[6]),
                    fun.components[3](centr + hlgth * XGK[6])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[8]),
                    fun.components[1](centr + hlgth * XGK[8]), fun.components[2](centr + hlgth * XGK[8]),
                    fun.components[3](centr + hlgth * XGK[8])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[10]),
                    fun.components[1](centr + hlgth * XGK[10]), fun.components[2](centr + hlgth * XGK[10]),
                    fun.components[3](centr + hlgth * XGK[10])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[12]),
                    fun.components[1](centr + hlgth * XGK[12]), fun.components[2](centr + hlgth * XGK[12]),
                    fun.components[3](centr + hlgth * XGK[12])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[14]),
                    fun.components[1](centr + hlgth * XGK[14]), fun.components[2](centr + hlgth * XGK[14]),
                    fun.components[3](centr + hlgth * XGK[14])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[16]),
                    fun.components[1](centr + hlgth * XGK[16]), fun.components[2](centr + hlgth * XGK[16]),
                    fun.components[3](centr + hlgth * XGK[16])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[18]),
                    fun.components[1](centr + hlgth * XGK[18]), fun.components[2](centr + hlgth * XGK[18]),
                    fun.components[3](centr + hlgth * XGK[18])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[20]),
                    fun.components[1](centr + hlgth * XGK[20]), fun.components[2](centr + hlgth * XGK[20]),
                    fun.components[3](centr + hlgth * XGK[20])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[22]),
                    fun.components[1](centr + hlgth * XGK[22]), fun.components[2](centr + hlgth * XGK[22]),
                    fun.components[3](centr + hlgth * XGK[22])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[24]),
                    fun.components[1](centr + hlgth * XGK[24]), fun.components[2](centr + hlgth * XGK[24]),
                    fun.components[3](centr + hlgth * XGK[24])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[26]),
                    fun.components[1](centr + hlgth * XGK[26]), fun.components[2](centr + hlgth * XGK[26]),
                    fun.components[3](centr + hlgth * XGK[26])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[28]),
                    fun.components[1](centr + hlgth * XGK[28]), fun.components[2](centr + hlgth * XGK[28]),
                    fun.components[3](centr + hlgth * XGK[28])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[30]),
                    fun.components[1](centr + hlgth * XGK[30]), fun.components[2](centr + hlgth * XGK[30]),
                    fun.components[3](centr + hlgth * XGK[30])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[32]),
                    fun.components[1](centr + hlgth * XGK[32]), fun.components[2](centr + hlgth * XGK[32]),
                    fun.components[3](centr + hlgth * XGK[32])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[34]),
                    fun.components[1](centr + hlgth * XGK[34]), fun.components[2](centr + hlgth * XGK[34]),
                    fun.components[3](centr + hlgth * XGK[34])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[36]),
                    fun.components[1](centr + hlgth * XGK[36]), fun.components[2](centr + hlgth * XGK[36]),
                    fun.components[3](centr + hlgth * XGK[36])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[38]),
                    fun.components[1](centr + hlgth * XGK[38]), fun.components[2](centr + hlgth * XGK[38]),
                    fun.components[3](centr + hlgth * XGK[38])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[40]),
                    fun.components[1](centr + hlgth * XGK[40]), fun.components[2](centr + hlgth * XGK[40]),
                    fun.components[3](centr + hlgth * XGK[40])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[42]),
                    fun.components[1](centr + hlgth * XGK[42]), fun.components[2](centr + hlgth * XGK[42]),
                    fun.components[3](centr + hlgth * XGK[42])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[44]),
                    fun.components[1](centr + hlgth * XGK[44]), fun.components[2](centr + hlgth * XGK[44]),
                    fun.components[3](centr + hlgth * XGK[44])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[46]),
                    fun.components[1](centr + hlgth * XGK[46]), fun.components[2](centr + hlgth * XGK[46]),
                    fun.components[3](centr + hlgth * XGK[46])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[48]),
                    fun.components[1](centr + hlgth * XGK[48]), fun.components[2](centr + hlgth * XGK[48]),
                    fun.components[3](centr + hlgth * XGK[48])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[50]),
                    fun.components[1](centr + hlgth * XGK[50]), fun.components[2](centr + hlgth * XGK[50]),
                    fun.components[3](centr + hlgth * XGK[50])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[52]),
                    fun.components[1](centr + hlgth * XGK[52]), fun.components[2](centr + hlgth * XGK[52]),
                    fun.components[3](centr + hlgth * XGK[52])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[54]),
                    fun.components[1](centr + hlgth * XGK[54]), fun.components[2](centr + hlgth * XGK[54]),
                    fun.components[3](centr + hlgth * XGK[54])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[56]),
                    fun.components[1](centr + hlgth * XGK[56]), fun.components[2](centr + hlgth * XGK[56]),
                    fun.components[3](centr + hlgth * XGK[56])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[58]),
                    fun.components[1](centr + hlgth * XGK[58]), fun.components[2](centr + hlgth * XGK[58]),
                    fun.components[3](centr + hlgth * XGK[58])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[60]),
                    fun.components[1](centr + hlgth * XGK[60]), fun.components[2](centr + hlgth * XGK[60]),
                    fun.components[3](centr + hlgth * XGK[60])])],




     [f64x4::from_array([fun.components[0](centr + hlgth * XGK[1]),
                 fun.components[1](centr + hlgth * XGK[1]), fun.components[2](centr + hlgth * XGK[1]),
                 fun.components[3](centr + hlgth * XGK[1])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[3]),
                    fun.components[1](centr + hlgth * XGK[3]), fun.components[2](centr + hlgth * XGK[3]),
                    fun.components[3](centr + hlgth * XGK[3])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[5]),
                    fun.components[1](centr + hlgth * XGK[5]), fun.components[2](centr + hlgth * XGK[5]),
                    fun.components[3](centr + hlgth * XGK[5])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[7]),
                    fun.components[1](centr + hlgth * XGK[7]), fun.components[2](centr + hlgth * XGK[7]),
                    fun.components[3](centr + hlgth * XGK[7])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[9]),
                    fun.components[1](centr + hlgth * XGK[9]), fun.components[2](centr + hlgth * XGK[9]),
                    fun.components[3](centr + hlgth * XGK[9])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[11]),
                    fun.components[1](centr + hlgth * XGK[11]), fun.components[2](centr + hlgth * XGK[11]),
                    fun.components[3](centr + hlgth * XGK[11])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[13]),
                    fun.components[1](centr + hlgth * XGK[13]), fun.components[2](centr + hlgth * XGK[13]),
                    fun.components[3](centr + hlgth * XGK[13])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[15]),
                    fun.components[1](centr + hlgth * XGK[15]), fun.components[2](centr + hlgth * XGK[15]),
                    fun.components[3](centr + hlgth * XGK[15])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[17]),
                    fun.components[1](centr + hlgth * XGK[17]), fun.components[2](centr + hlgth * XGK[17]),
                    fun.components[3](centr + hlgth * XGK[17])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[19]),
                    fun.components[1](centr + hlgth * XGK[19]), fun.components[2](centr + hlgth * XGK[19]),
                    fun.components[3](centr + hlgth * XGK[19])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[21]),
                    fun.components[1](centr + hlgth * XGK[21]), fun.components[2](centr + hlgth * XGK[21]),
                    fun.components[3](centr + hlgth * XGK[21])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[23]),
                    fun.components[1](centr + hlgth * XGK[23]), fun.components[2](centr + hlgth * XGK[23]),
                    fun.components[3](centr + hlgth * XGK[23])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[25]),
                    fun.components[1](centr + hlgth * XGK[25]), fun.components[2](centr + hlgth * XGK[25]),
                    fun.components[3](centr + hlgth * XGK[25])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[27]),
                    fun.components[1](centr + hlgth * XGK[27]), fun.components[2](centr + hlgth * XGK[27]),
                    fun.components[3](centr + hlgth * XGK[27])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[29]),
                    fun.components[1](centr + hlgth * XGK[29]), fun.components[2](centr + hlgth * XGK[29]),
                    fun.components[3](centr + hlgth * XGK[29])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[31]),
                    fun.components[1](centr + hlgth * XGK[31]), fun.components[2](centr + hlgth * XGK[31]),
                    fun.components[3](centr + hlgth * XGK[31])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[33]),
                    fun.components[1](centr + hlgth * XGK[33]), fun.components[2](centr + hlgth * XGK[33]),
                    fun.components[3](centr + hlgth * XGK[33])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[35]),
                    fun.components[1](centr + hlgth * XGK[35]), fun.components[2](centr + hlgth * XGK[35]),
                    fun.components[3](centr + hlgth * XGK[35])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[37]),
                    fun.components[1](centr + hlgth * XGK[37]), fun.components[2](centr + hlgth * XGK[37]),
                    fun.components[3](centr + hlgth * XGK[37])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[39]),
                    fun.components[1](centr + hlgth * XGK[39]), fun.components[2](centr + hlgth * XGK[39]),
                    fun.components[3](centr + hlgth * XGK[39])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[41]),
                    fun.components[1](centr + hlgth * XGK[41]), fun.components[2](centr + hlgth * XGK[41]),
                    fun.components[3](centr + hlgth * XGK[41])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[43]),
                    fun.components[1](centr + hlgth * XGK[43]), fun.components[2](centr + hlgth * XGK[43]),
                    fun.components[3](centr + hlgth * XGK[43])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[45]),
                    fun.components[1](centr + hlgth * XGK[45]), fun.components[2](centr + hlgth * XGK[45]),
                    fun.components[3](centr + hlgth * XGK[45])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[47]),
                    fun.components[1](centr + hlgth * XGK[47]), fun.components[2](centr + hlgth * XGK[47]),
                    fun.components[3](centr + hlgth * XGK[47])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[49]),
                    fun.components[1](centr + hlgth * XGK[49]), fun.components[2](centr + hlgth * XGK[49]),
                    fun.components[3](centr + hlgth * XGK[49])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[51]),
                    fun.components[1](centr + hlgth * XGK[51]), fun.components[2](centr + hlgth * XGK[51]),
                    fun.components[3](centr + hlgth * XGK[51])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[53]),
                    fun.components[1](centr + hlgth * XGK[53]), fun.components[2](centr + hlgth * XGK[53]),
                    fun.components[3](centr + hlgth * XGK[53])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[55]),
                    fun.components[1](centr + hlgth * XGK[55]), fun.components[2](centr + hlgth * XGK[55]),
                    fun.components[3](centr + hlgth * XGK[55])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[57]),
                    fun.components[1](centr + hlgth * XGK[57]), fun.components[2](centr + hlgth * XGK[57]),
                    fun.components[3](centr + hlgth * XGK[57])]),
         f64x4::from_array([fun.components[0](centr + hlgth * XGK[59]),
                    fun.components[1](centr + hlgth * XGK[59]), fun.components[2](centr + hlgth * XGK[59]),
                    fun.components[3](centr + hlgth * XGK[59])]),
         f64x4::from_array([0.0,0.0,0.0,0.0])])


}


/*
pub fn fvec_simd2(fun: &FnVec4, centr : f64, hlgth : f64) -> Vec<f64x4>{
    let mut vec = vec![];
    for x in XGK {
        vec.push(f64x4::new(fun.components[0](centr + hlgth * x),
            fun.components[1](centr + hlgth * x), fun.components[2](centr + hlgth * x),
            fun.components[3](centr + hlgth * x)));
    }
    vec.push(f64x4::splat(0.0));
    vec

}


 */
#[cfg(test)]
mod tests {
    use std::simd::{f64x4, Simd};
    use std::time::Instant;
    use crate::funct_vector::FnVec4;
    use crate::qk61::Qk61;
    use crate::qk61_1dvec3::Qk611DVec3;
    use crate::qk61_4vec_simd2::Qk61Vec4Simd2;
    use crate::qk61_4vec_simd3::Qk61Vec4Simd3;
    use crate::qk61_4vec_simd4::Qk61Vec4Simd4;
    use crate::qk61_simd::Qk61Simd;
    use crate::qk61_4vec_simd::Qk61Vec4Simd;
    use crate::qk61_simd2::Qk61Simd2;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f = |x:f64| x.cos();

        let a = 0.0;
        let b = 1.0;
        let qks = Qk61Vec4Simd3{};
        let qks2 = Qk61Vec4Simd4{};
        let qk = Qk61Simd {};
        let qk2 = Qk61Simd2 {};
        let fun = FnVec4{ components : [Box::new(f),Box::new(f),Box::new(f),Box::new(f)]};

        let null = f64x4::splat(0.0);
        let null2 = packed_simd_2::f64x4::splat(0.0);
        //let mut res_simd = (null,null,null,null);
        let mut res_simd = (null2,null2,null2,null2);
        let mut res_simd2 = (null,null,null,null);

        for k in 0..100 {
            let start = Instant::now();
            let res1 = qk.integrate(&f, a, b);
            let res2 = qk.integrate(&f,a,b);
            let res3 = qk.integrate(&f,a,b);
            let res4 = qk.integrate(&f,a,b);
            println!("normal {:?}", start.elapsed());
            let start = Instant::now();
            let res1 = qk2.integrate(&f, a, b);
            let res2 = qk2.integrate(&f,a,b);
            let res3 = qk2.integrate(&f,a,b);
            let res4 = qk2.integrate(&f,a,b);
            let rres = (Simd::from_array([res1.0,res2.0,res3.0,res4.0]),
                        Simd::from_array([res1.1,res2.1,res3.1,res4.1]),
                        Simd::from_array([res1.2,res2.2,res3.2,res4.2]),
                        Simd::from_array([res1.3,res2.3,res3.3,res4.3]));
            println!("normal2 {:?}", start.elapsed());
            let start = Instant::now();
            res_simd = qks.integrate(&fun, a, b);
            println!("simd {:?}", start.elapsed());
            let start = Instant::now();
            res_simd2 = qks2.integrate(&fun, a, b);
            println!("simd2 {:?}", start.elapsed());
        }
        println!("simd {:?}",res_simd);
        println!("simd2 {:?}",res_simd2);
    }
}




