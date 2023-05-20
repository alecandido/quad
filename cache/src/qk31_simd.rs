use std::simd::{Simd, SimdFloat};
use crate::qk::*;


pub struct Qk31Simd {}
///     Parameters:
///
///     On entry:
///         f   :   f64
///                 function
///
///         a   :   f64
///                 lower limit of integration
///
///         b   :   f64
///                 upper limit of integration
///
///     On return:
///         result  :   f64
///                     approximation to the integral i result is computed by applying
///                     the 15-point kronrod rule (resk) obtained by optimal addition
///                     of abscissae to the7-point gauss rule(resg).
///
///         abserr  :   f64
///                     estimate of the modulus of the absolute error, which should not
///                     exceed abs(i-result)
///
///         resabs  :   f64
///                     approximation to the integral j
///
///         resasc  :   f64
///                     approximation to the integral of abs(f-i/(b-a)) over (a,b)
///
///     The abscissae and weights are given for the interval (-1,1).
///     Because of symmetry only the positive abscissae and their
///     corresponding weights are given.
///
///         xgk     :   abscissae of the 15-point kronrod rule
///                     xgk(2), xgk(4), ...  abscissae of the 7-point
///                     gauss rule
///                     xgk(1), xgk(3), ...  abscissae which are optimally
///                     added to the 7-point gauss rule
///
///         wgk     :   weights of the 15-point kronrod rule
///
///         wg      :   weights of the 7-point gauss rule
///
///
///     Gauss quadrature weights and kronrod quadrature abscissae and weights
///     as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
///     bell labs, nov. 1981.
///
///
///


const XGK: Simd<f64,32> = Simd::from_array([-0.998002298693397060285172840152271, -0.987992518020485428489565718586613,
    -0.967739075679139134257347978784337, -0.937273392400705904307758947710209,
    -0.897264532344081900882509656454496, -0.848206583410427216200648320774217,
    -0.790418501442465932967649294817947, -0.724417731360170047416186054613938,
    -0.650996741297416970533735895313275, -0.570972172608538847537226737253911,
    -0.485081863640239680693655740232351, -0.394151347077563369897207370981045,
    -0.299180007153168812166780024266389, -0.201194093997434522300628303394596,
    -0.101142066918717499027074231447392, 0.000000000000000000000000000000000,
    0.101142066918717499027074231447392, 0.201194093997434522300628303394596,
    0.299180007153168812166780024266389, 0.394151347077563369897207370981045,
    0.485081863640239680693655740232351, 0.570972172608538847537226737253911,
    0.650996741297416970533735895313275, 0.724417731360170047416186054613938,
    0.790418501442465932967649294817947, 0.848206583410427216200648320774217,
    0.897264532344081900882509656454496, 0.937273392400705904307758947710209,
    0.967739075679139134257347978784337, 0.987992518020485428489565718586613,
    0.998002298693397060285172840152271, 0.000000000000000000000000000000000]);


const WGK: Simd<f64,32> = Simd::from_array([0.005377479872923348987792051430128, 0.015007947329316122538374763075807,
    0.025460847326715320186874001019653, 0.035346360791375846222037948478360,
    0.044589751324764876608227299373280, 0.053481524690928087265343147239430,
    0.062009567800670640285139230960803, 0.069854121318728258709520077099147,
    0.076849680757720378894432777482659, 0.083080502823133021038289247286104,
    0.088564443056211770647275443693774, 0.093126598170825321225486872747346,
    0.096642726983623678505179907627589, 0.099173598721791959332393173484603,
    0.100769845523875595044946662617570, 0.101330007014791549017374792767493,
    0.100769845523875595044946662617570, 0.099173598721791959332393173484603,
    0.096642726983623678505179907627589, 0.093126598170825321225486872747346,
    0.088564443056211770647275443693774, 0.083080502823133021038289247286104,
    0.076849680757720378894432777482659, 0.069854121318728258709520077099147,
    0.062009567800670640285139230960803, 0.053481524690928087265343147239430,
    0.044589751324764876608227299373280, 0.035346360791375846222037948478360,
    0.025460847326715320186874001019653, 0.015007947329316122538374763075807,
    0.005377479872923348987792051430128, 0.000000000000000000000000000000000]);



const WG: Simd<f64,32> = Simd::from_array([0.0, 0.03075324199611726835462839357720, 0.0, 0.07036604748810812470926741645066,
    0.0, 0.10715922046717193501186954668586, 0.0, 0.13957067792615431444780479451102,
    0.0, 0.16626920581699393355320086048120, 0.0, 0.18616100001556221102680056186642,
    0.0, 0.19843148532711157645611832644383, 0.0, 0.20257824192556127288062019996751,
    0.0, 0.19843148532711157645611832644383, 0.0, 0.18616100001556221102680056186642,
    0.0, 0.16626920581699393355320086048120, 0.0, 0.13957067792615431444780479451102,
    0.0, 0.10715922046717193501186954668586, 0.0, 0.07036604748810812470926741645066,
    0.0, 0.030753241996117268354628393577204, 0.0, 0.000000000000000000000000000000000]);






impl Qk for Qk31Simd {
    fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64) -> (f64, f64, f64, f64){
        let hlgth: f64 = 0.5 * (b - a);
        let dhlgth: f64 = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let fv = fvec_simd(f, centr, hlgth);
        let resk = (fv * WGK).reduce_sum();
        let reskh = resk * 0.5;
        let reskhs = Simd::from_array([reskh;32]);
        let mut resabs = (fv.abs() * WGK).reduce_sum();
        let resg = ( fv * WG).reduce_sum();
        let mut resasc = ( WGK * ( fv - reskhs).abs()).reduce_sum();

        let result = resk * hlgth;
        resabs = resabs * dhlgth;
        resasc = resasc * dhlgth;
        let mut abserr = ((resk - resg) * hlgth).abs();
        if resasc != 0.0 && abserr != 0.0 {
            abserr = resasc * 1.0_f64.min((200.0 * abserr / resasc).powf(1.5));
        }
        if resabs > UFLOW / (50.0 * EPMACH) {
            abserr = abserr.max((EPMACH * 50.0) * resabs);
        }

        (result, abserr, resabs, resasc)
    }
}



pub fn fvec_simd(f : &dyn Fn(f64)->f64, centr : f64, hlgth : f64) -> Simd<f64,32>{
    Simd::from_array([f(centr + hlgth * XGK[0]), f(centr + hlgth * XGK[1]), f(centr + hlgth * XGK[2]),
        f(centr + hlgth * XGK[3]), f(centr + hlgth * XGK[4]), f(centr + hlgth * XGK[5]),
        f(centr + hlgth * XGK[6]), f(centr + hlgth * XGK[7]), f(centr + hlgth * XGK[8]),
        f(centr + hlgth * XGK[9]), f(centr + hlgth * XGK[10]),f(centr + hlgth * XGK[11]),
        f(centr + hlgth * XGK[12]), f(centr + hlgth * XGK[13]), f(centr + hlgth * XGK[14]),
        f(centr + hlgth * XGK[15]), f(centr + hlgth * XGK[16]), f(centr + hlgth * XGK[17]),
        f(centr + hlgth * XGK[18]), f(centr + hlgth * XGK[19]), f(centr + hlgth * XGK[20]),
        f(centr + hlgth * XGK[21]), f(centr + hlgth * XGK[22]), f(centr + hlgth * XGK[23]),
        f(centr + hlgth * XGK[24]), f(centr + hlgth * XGK[25]),f(centr + hlgth * XGK[26]),
        f(centr + hlgth * XGK[27]), f(centr + hlgth * XGK[28]), f(centr + hlgth * XGK[29]),
        f(centr + hlgth * XGK[30]), 0.000000000000000000000000000000000])
}

#[cfg(test)]
mod tests {
    use std::time::Instant;
    use crate::qk15::Qk15;
    use crate::qk15_simd::Qk15Simd;
    use crate::qk21::Qk21;
    use crate::qk21_simd::Qk21Simd;
    use crate::qk31::Qk31;
    use crate::qk31_simd::Qk31Simd;
    use crate::qk51::Qk51;
    use crate::qk51_simd::Qk51Simd;
    use crate::qk61_simd::Qk61Simd;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f = |x:f64| x.cos();
        let a = 0.0;
        let b = 1.0;
        let qks = Qk31Simd {};
        let qk = Qk31{};

        let mut res2 = (0.0,0.0,0.0,0.0);
        let mut res1 = res2.clone();

        for k in 0..100 {
            let start = Instant::now();
            res2 = qk.integrate(&f, a, b);
            println!("normal {:?}", start.elapsed());
            let start = Instant::now();
            res1 = qks.integrate(&f, a, b);
            println!("simd {:?}", start.elapsed());
        }
        println!("{:?}",res1);
        println!("{:?}",res2);
    }
}







