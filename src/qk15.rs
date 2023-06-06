use crate::qk_array::qk_array_quadrature;
use crate::qk_vec::qk_quadrature_vec;

pub fn qk15_array_quadrature<const N:usize,F>(f: F, a: f64, b: f64) -> ([f64; N], f64, f64)
    where F : Fn(f64) -> [f64; N]
{
    qk_array_quadrature(f, a, b, &XGK15, &WGK15, &WG15)
}

//pub fn qk15_quadrature_v<F>(f: F, a: f64, b: f64) -> (Vec<f64>, f64, f64)
//    where F : Fn(f64) -> Vec<f64>
//{
//    qk_quadrature_v(f,a,b,&XGK15,&WGK15,&WG15)
//}
pub fn qk15_vec_quadrature<F>(f: F, a: f64, b: f64) -> (Vec<f64>, f64, f64)
    where F : Fn(f64) -> Vec<f64>
{
    qk_quadrature_vec(f, a, b, &XGK15, &WGK15, &WG15)
}


const XGK15 : [f64;7] = [0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144838258730, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245];

const WGK15 : [f64;8] = [0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714];

const WG15 : [f64;4] = [0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
    0.381830050505118944950369775488975, 0.417959183673469387755102040816327];



#[cfg(test)]
mod tests {
    use std::time::Instant;
    use crate::qk15::{qk15_array_quadrature, qk15_vec_quadrature};

    #[test]
    fn test(){
        let a = 0.0;
        let b = 1.0;
        let f_array = |x:f64| [x.cos(),x.sin()];
        let f_vec = |x:f64| vec![x.cos(),x.sin()];

        let max = 100;

        let mut res_array;
        let mut res_vec;

        let (mut t1,mut t2) = (0.0,0.0);

        for i in 0..max {
            let time = Instant::now();
            res_array = qk15_array_quadrature(&f_array,a,b);
            if i > 10 { t1 += time.elapsed().as_secs_f64(); }

            let time = Instant::now();
            res_vec = qk15_vec_quadrature(&f_vec,a,b);
            if i > 10 { t2 += time.elapsed().as_secs_f64(); }

            if i == max -1 {
                println!("array res : {:?}", res_array);
                println!("vec res : {:?}", res_vec);
            }
        }

        t1 /= max as f64 - 10.0;
        t2 /= max as f64 - 10.0;

        println!("array time : {t1}; vec time : {t2}");

    }
}
