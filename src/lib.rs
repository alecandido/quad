pub mod constants;
pub mod qag_par_integrator_result;
pub mod result_state;
pub mod qk_array;
pub mod qk61;
pub mod qk51;
pub mod qk41;
pub mod qk31;
pub mod qk21;
pub mod qk15;
pub mod qag_par_integration_result;
pub mod semi_infinite_function;
pub mod qag_array_integration_result;
pub mod qag_array_integrator_result;
pub mod qag_array;
mod qk_vec;
mod qag_vec_integration_result;
mod qag_vec_integrator_result;
mod qag_vec;
mod qag_par_iter;

use pyo3::prelude::*;
use crate::constants::Myf64;
use crate::qag_array::QagArray;
use crate::qag_array_integration_result::QagArrayIntegrationResult;
use crate::qag_vec::QagVec;

fn lambda_eval<const N: usize>(ob : &PyAny, z : f64) -> [f64; N]{
    let f = |x:f64| {
        let y = (x,);
        ob.call1(y).unwrap().extract::<[f64; N]>().unwrap()
    };
    f(z)
}

fn lambda_eval_vec(ob : &PyAny, z : f64) -> Vec<f64>{
    let f = |x:f64| {
        let y = (x,);
        ob.call1(y).unwrap().extract::<Vec<f64>>().unwrap()
    };
    f(z)
}





#[pyfunction]
fn qag_array(ob : &PyAny, a: f64, b: f64, epsabs : Option<f64>, epsrel : Option<f64>, key : Option<i32>, limit : Option<usize>, dim : usize,
             points : Option<Vec<f64>>, more_info : Option<bool>)
             -> Py<QagsResult> {

    let pointss = {
        if points.is_some() {
            points.unwrap()
        }
        else {
            [0.0;0].to_vec()
        }
    };
    let limitt = {
        if limit.is_some(){ limit.unwrap() }
        else { 1000 }
    };
    let keyy = {
        if key.is_some(){ key.unwrap() }
        else{ 2 }
    };
    let epsabss = {
        if epsabs.is_some(){ epsabs.unwrap() }
        else { 1e-4 }
    };
    let epsrell = {
        if epsrel.is_some(){ epsrel.unwrap() }
        else { 1e-4 }
    };
    let more_infoo = {
        if more_info.is_some() { more_info.unwrap() }
        else { false }
    };



    let out = Python::with_gil(|py| {
        let qag = QagArray {key : keyy,limit : limitt, points : pointss, more_info : more_infoo};

        match dim {
            1 => {
                let f = |x: f64| lambda_eval::<1>(ob, x);
                let (result, abserr, more_info) = qag_extract(&f, &qag, a, b, epsabss, epsrell);
                Py::new(py, QagsResult {result,abserr,more_info}).unwrap()

            },
            2 => {
                let f = |x: f64| lambda_eval::<2>(ob, x);
                let (result, abserr, more_info) = qag_extract(&f, &qag, a, b, epsabss, epsrell);
                Py::new(py, QagsResult {result,abserr,more_info}).unwrap()
            },

            3 => {
                let f = |x: f64| lambda_eval::<3>(ob, x);
                let (result, abserr, more_info) = qag_extract(&f, &qag, a, b, epsabss, epsrell);
                Py::new(py, QagsResult {result,abserr,more_info}).unwrap()
            },

            4 => {
                let f = |x: f64| lambda_eval::<4>(ob, x);
                let (result, abserr, more_info) = qag_extract(&f, &qag, a, b, epsabss, epsrell);
                Py::new(py, QagsResult {result,abserr,more_info}).unwrap()
            },

            5 => {
                let f = |x: f64| lambda_eval::<5>(ob, x);
                let (result, abserr, more_info) = qag_extract(&f, &qag, a, b, epsabss, epsrell);
                Py::new(py, QagsResult {result,abserr,more_info}).unwrap()
            },

            6 => {
                let f = |x: f64| lambda_eval::<6>(ob, x);
                let (result, abserr, more_info) = qag_extract(&f, &qag, a, b, epsabss, epsrell);
                Py::new(py, QagsResult {result,abserr,more_info}).unwrap()
            },

            7 => {
                let f = |x: f64| lambda_eval::<7>(ob, x);
                let (result, abserr, more_info) = qag_extract(&f, &qag, a, b, epsabss, epsrell);
                Py::new(py, QagsResult {result,abserr,more_info}).unwrap()
            },

            8 => {
                let f = |x: f64| lambda_eval::<8>(ob, x);
                let (result, abserr, more_info) = qag_extract(&f, &qag, a, b, epsabss, epsrell);
                Py::new(py, QagsResult {result,abserr,more_info}).unwrap()
            },


            _ => {
                let f = |x: f64| lambda_eval::<1>(ob, x);
                let (result, abserr, more_info) = qag_extract(&f, &qag, a, b, epsabss, epsrell);
                Py::new(py, QagsResult {result,abserr,more_info}).unwrap()
            },
        }
    });
    out
}

#[pyfunction]
fn qag_vec(ob : &PyAny, a: f64, b: f64, epsabs : Option<f64>, epsrel : Option<f64>, key : Option<i32>, limit : Option<usize>,
             points : Option<Vec<f64>>, more_info : Option<bool>)
             -> Py<QagsResult> {

    let pointss = {
        if points.is_some() {
            points.unwrap()
        }
        else {
            [0.0;0].to_vec()
        }
    };
    let limitt = {
        if limit.is_some(){ limit.unwrap() }
        else { 1000 }
    };
    let keyy = {
        if key.is_some(){ key.unwrap() }
        else{ 2 }
    };
    let epsabss = {
        if epsabs.is_some(){ epsabs.unwrap() }
        else { 1e-4 }
    };
    let epsrell = {
        if epsrel.is_some(){ epsrel.unwrap() }
        else { 1e-4 }
    };
    let more_infoo = {
        if more_info.is_some() { more_info.unwrap() }
        else { false }
    };



    let out = Python::with_gil(|py| {
        let qag = QagVec{key : keyy,limit : limitt, points : pointss, more_info : more_infoo};

        let f = |x: f64| lambda_eval_vec(ob, x);
        let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
        let (result,abserr,more_inf) = (res.result,res.abserr,res.more_info);
        if more_inf.is_none(){
            Py::new(py, QagsResult {result,abserr,more_info : None}).unwrap()
        }
        else {
            let mut more_inf_py : Vec<(f64,f64,f64,Vec<f64>)> = vec![];
            let more_inf_unwrapped = more_inf.unwrap();
            let (mut hash, mut heap) = (more_inf_unwrapped.hash, more_inf_unwrapped.heap);
            let (neval, last) = (more_inf_unwrapped.neval, more_inf_unwrapped.last);
            for _k in 0..heap.len(){
                let old_interval = heap.pop().unwrap();
                let ((x,y),old_err) = (old_interval.interval,old_interval.err);
                let old_res = hash.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
                more_inf_py.push((x,y,old_err,old_res));
            }
            Py::new(py, QagsResult {result,abserr,more_info : Some((neval,last,more_inf_py))}).unwrap()

        }
    });
    out
}

pub fn qag_extract<const N : usize,F>(f : &F, qag : &QagArray, a : f64, b : f64, epsabss : f64, epsrell : f64)
                                      -> (Vec<f64>,f64,Option<(i32,usize,Vec<(f64,f64,f64,Vec<f64>)>)>)
    where F : Fn(f64) -> [f64;N]{

    let res : QagArrayIntegrationResult<N> = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
    let (result,abserr,more_inf) = (res.result.to_vec(),res.abserr,res.more_info);
    if more_inf.is_none(){
        (result,abserr,None)
    }
    else {
        let mut more_inf_py : Vec<(f64,f64,f64,Vec<f64>)> = vec![];
        let more_inf_unwrapped = more_inf.unwrap();
        let (mut hash, mut heap) = (more_inf_unwrapped.hash, more_inf_unwrapped.heap);
        let (neval, last) = (more_inf_unwrapped.neval, more_inf_unwrapped.last);
        for _k in 0..heap.len(){
            let old_interval = heap.pop().unwrap();
            let ((x,y),old_err) = (old_interval.interval,old_interval.err);
            let old_res = hash.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
            more_inf_py.push((x,y,old_err,old_res.to_vec()));
        }
        (result,abserr,Some((neval,last,more_inf_py)))
    }
}

#[pyclass]
struct QagsResult {
    #[pyo3(get, set)]
    pub result : Vec<f64>,
    #[pyo3(get, set)]
    pub abserr : f64,
    #[pyo3(get, set)]
    pub more_info : Option<(i32,usize,Vec<(f64,f64,f64,Vec<f64>)>)>
}

#[pyclass]
struct MoreInfoPy {
    neval : i32,
    last : usize,
    hash : Vec<(f64,f64,f64,Vec<f64>)>,
}




/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn quad(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(qag_array, m)?)?;
    m.add_function(wrap_pyfunction!(qag_vec, m)?)?;
    Ok(())
}






#[cfg(test)]
mod tests {
    use std::time::Instant;

    #[test]
    fn array_vs_vec(){

        for _i in 0..100{
            let start = Instant::now();
            const MAX: usize = 128;
            let mut array = [0.0; MAX];
            for k in 0..MAX {
                array[k] = k as f64;
            }
            println!("array time : {:?}",start.elapsed());

            let start = Instant::now();
            let mut vec = Vec::with_capacity(MAX);
            for k in 0..MAX {
                vec.push(k as f64);
            }
            println!("vec time : {:?}",start.elapsed());

            let start = Instant::now();
            let mut vec = vec![0.0; MAX];
            for k in 0..MAX {
                vec[k] = k as f64;
            }
            println!("vec init time : {:?}",start.elapsed());

            println!("{:?}",array.iter().sum::<f64>());
            println!("{:?}",vec.iter().sum::<f64>());
        }



    }




}



