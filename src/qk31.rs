use crate::qk::*;

pub struct Qk31{}

const XGK : [f64;16] = [0.998002298693397060285172840152271, 0.987992518020485428489565718586613,
                        0.967739075679139134257347978784337, 0.937273392400705904307758947710209,
                        0.897264532344081900882509656454496, 0.848206583410427216200648320774217,
                        0.790418501442465932967649294817947, 0.724417731360170047416186054613938,
                        0.650996741297416970533735895313275, 0.570972172608538847537226737253911,
                        0.485081863640239680693655740232351, 0.394151347077563369897207370981045,
                        0.299180007153168812166780024266389, 0.201194093997434522300628303394596,
                        0.101142066918717499027074231447392, 0.000000000000000000000000000000000];


const WGK : [f64;16] = [0.005377479872923348987792051430128, 0.015007947329316122538374763075807,
                        0.025460847326715320186874001019653, 0.035346360791375846222037948478360,
                        0.044589751324764876608227299373280, 0.053481524690928087265343147239430,
                        0.062009567800670640285139230960803, 0.069854121318728258709520077099147,
                        0.076849680757720378894432777482659, 0.083080502823133021038289247286104,
                        0.088564443056211770647275443693774, 0.093126598170825321225486872747346,
                        0.096642726983623678505179907627589, 0.099173598721791959332393173484603,
                        0.100769845523875595044946662617570, 0.101330007014791549017374792767493];

const WG : [f64;8] = [0.030753241996117268354628393577204, 0.070366047488108124709267416450667,
                      0.107159220467171935011869546685869, 0.139570677926154314447804794511028,
                      0.166269205816993933553200860481209, 0.186161000015562211026800561866423,
                      0.198431485327111576456118326443839, 0.202578241925561272880620199967519];

impl Qk for Qk31 {
    fn integrate(&self, f: &dyn Fn(f64) -> f64, a: f64, b: f64, ) -> (f64, f64, f64, f64) {
        let hlgth: f64 = 0.5 * (b - a);
        let dhlgth: f64 = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let mut fv1: Vec<f64> = vec![0.0; 15];
        let mut fv2: Vec<f64> = vec![0.0; 15];

        //compute the 31-point kronrod approximation to
        //the integral, and estimate the absolute error.

        let fc : f64 = f(centr);
        let mut resg = WG[7] * fc;
        let mut resk = WGK[15] * fc;
        let mut resabs = resk.abs();


        for j in 1..8 {
            let jtw = 2 * j;
            let absc = hlgth * XGK[jtw - 1];
            let fval1 : f64 = f(centr - absc);
            let fval2 : f64 = f(centr + absc);
            fv1[jtw - 1] = fval1;
            fv2[jtw - 1] = fval2;
            let fsum = fval1 + fval2;
            resg += WG[j - 1] * fsum;
            resk += WGK[jtw - 1] * fsum;
            resabs += WGK[jtw - 1] * (fval1.abs() + fval2.abs());
        }


        for j in 1..9 {
            let jtwm1 = 2 * j - 1;
            let absc = hlgth * XGK[jtwm1 - 1];
            let fval1 : f64 = f(centr - absc);
            let fval2 : f64 = f(centr + absc);
            fv1[jtwm1 - 1] = fval1;
            fv2[jtwm1 - 1] = fval2;
            let fsum = fval1 + fval2;
            resk += WGK[jtwm1 - 1] * fsum;
            resabs += WGK[jtwm1 - 1] * (fval1.abs() + fval2.abs());
        }

        let reskh = resk * 0.5;
        let mut resasc = WGK[15] * (fc - reskh).abs();

        for j in 1..16 {
            resasc += WGK[j - 1] * ((fv1[j - 1] - reskh).abs() + (fv2[j - 1] - reskh).abs());
        }

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