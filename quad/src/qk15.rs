use crate::qk::qk_quadrature;
use crate::qk_ndarray::qk_quadrature_ndarray;

pub fn qk15_quadrature<F>(f: F, a: f64, b: f64) -> (Vec<f64>, f64, f64)
where
    F: Fn(f64) -> Vec<f64>,
{
    qk_quadrature(f, a, b, &XGK15, &WGK15, &WG15)
}

pub fn qk15_quadrature_ndarray<F>(f: F, a: f64, b: f64) -> (ndarray::Array1<f64>, f64, f64)
where
    F: Fn(f64) -> ndarray::Array1<f64>,
{
    qk_quadrature_ndarray(f, a, b, &XGK15, &WGK15, &WG15)
}

const XGK15: [f64; 7] = [
    0.991455371120812639206854697526329,
    0.949107912342758524526189684047851,
    0.864864423359769072789712788640926,
    0.741531185599394439863864773280788,
    0.586087235467691130294144838258730,
    0.405845151377397166906606412076961,
    0.207784955007898467600689403773245,
];

const WGK15: [f64; 8] = [
    0.022935322010529224963732008058970,
    0.063092092629978553290700663189204,
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714,
];

const WG15: [f64; 4] = [
    0.129484966168869693270611432679082,
    0.279705391489276667901467771423780,
    0.381830050505118944950369775488975,
    0.417959183673469387755102040816327,
];

#[cfg(test)]
mod tests {
    use crate::qk15::{qk15_quadrature, qk15_quadrature_ndarray};
    use ndarray::array;
    use std::time::Instant;

    #[test]
    fn vec_vs_ndarray() {
        let a = 0.0;
        let b = 10000.0;

        let f_vec = |x: f64| vec![x.sin(), x.cos()];
        let f_nd = |x: f64| array![x.sin(), x.cos()];

        for k in 0..25000 {
            let start = Instant::now();
            let res_vec = qk15_quadrature(&f_vec, a, b);
            let time_vec = start.elapsed().as_secs_f64();
            let start = Instant::now();
            let res_nd = qk15_quadrature_ndarray(&f_nd, a, b);
            let time_nd = start.elapsed().as_secs_f64();

            if k == 24999 {
                println!("{:?}", res_vec);
                println!("{:?}", res_nd);
            }

            println!("{:?}", time_vec);
            println!("{:?}", time_nd);
        }
    }
}
