mod constants;
mod qag_vec_norm_integration_result;
mod qag_vec_norm_integrator_result;
mod qage_vec_norm2;
mod qage_vec_norm3;
mod qage_vec_norm_parall;
mod qk61_vec_norm2;
mod result_state;


#[cfg(test)]
mod tests {
    use crate::qage_vec_norm3::QagVecNorm3;

    #[test]
    fn qag() {
        let a = 0.0;
        let b = 1.0e2;
        let key = 6;
        let limit = 1000000;
        let epsrel = 0.0;
        let epsabs = 15.3;

        let f = |x:f64| [x.cos(),x.sin(),x.sqrt(),x.sqrt()];
        let qag = QagVecNorm3{key,limit};

        let res = qag.qintegrate(&f,a,b,epsabs,epsrel).unwrap();
        println!("{:?}",res);




    }







}




