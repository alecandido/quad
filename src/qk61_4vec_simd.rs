use std::iter::zip;
use std::simd::Simd;
use std::time::Instant;
use crate::funct_vector::FnVec4;
use crate::qk::*;
use crate::qk61_1dvec_simd::Qk611DVec_Simd;

pub struct Qk61Vec4Simd {}

impl Qk61Vec4Simd {
    pub(crate) fn integrate(&self, fun: &FnVec4, a: f64, b: f64, ) -> (Simd<f64,4>, Simd<f64,4>, Simd<f64,4>, Simd<f64,4>,) {

        let qk61_1d = Qk611DVec_Simd{};
        let (mut result_vec,mut abserr_vec,mut resabs_vec,mut resasc_vec) =
            ([0.0;4],[0.0;4],[0.0;4],[0.0;4]);
        for k in 0..4 {
            (result_vec[k], abserr_vec[k], resabs_vec[k], resasc_vec[k]) = qk61_1d.integrate(&fun.components[k],a,b);
        }
        (Simd::from_array(result_vec),Simd::from_array(abserr_vec),
         Simd::from_array(resabs_vec),Simd::from_array(resasc_vec))
    }
}

#[cfg(test)]
mod tests {
    use std::time::Instant;
    use crate::funct_vector::FnVec4;
    use crate::qk61::Qk61;
    use crate::qk61_1dvec3::Qk611DVec3;
    use crate::qk61_1dvec_simd::Qk611DVec_Simd;
    use crate::qk61_4vec_simd::Qk61Vec4Simd;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f = |x:f64| x.cos();

        let a = 0.0;
        let b = 1.0;
        let qks = Qk61Vec4Simd{};
        let qk = Qk611DVec_Simd{};
        let fun = FnVec4{ components : [Box::new(f),Box::new(f),Box::new(f),Box::new(f)]};

        for k in 0..100 {
            let start = Instant::now();
            let res1 = qk.integrate(&f, a, b);
            let res2 = qk.integrate(&f,a,b);
            let res3 = qk.integrate(&f,a,b);
            let res4 = qk.integrate(&f,a,b);
            println!("normal {:?}", start.elapsed());
            let start = Instant::now();
            let res_simd = qks.integrate(&fun, a, b);
            println!("simd {:?}", start.elapsed());
        }
        //println!("{:?}",res1);
        //println!("{:?}",res1);
    }
}

