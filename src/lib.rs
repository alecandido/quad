pub mod constants;
pub mod result_state;
pub mod qk61;
pub mod qk51;
pub mod qk41;
pub mod qk31;
pub mod qk21;
pub mod qk15;
pub mod semi_infinite_function;
mod qk;
mod qag_integration_result;
mod qag_integrator_result;
mod qag;
mod qag_par;
use pyo3::prelude::*;
use crate::constants::Myf64;
use crate::qag::Qag;


fn lambda_eval(ob : &PyAny, z : f64) -> Vec<f64>{
    let f = |x:f64| {
        let y = (x,);
        ob.call1(y).unwrap().extract::<Vec<f64>>().unwrap()
    };
    f(z)
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
        let qag = Qag {key : keyy,limit : limitt, points : pointss, more_info : more_infoo};

        let f = |x: f64| lambda_eval(ob, x);
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


#[pyclass]
struct QagsResult {
    #[pyo3(get, set)]
    pub result : Vec<f64>,
    #[pyo3(get, set)]
    pub abserr : f64,
    #[pyo3(get, set)]
    pub more_info : Option<(i32,usize,Vec<(f64,f64,f64,Vec<f64>)>)>
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn quad(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(qag_vec, m)?)?;
    Ok(())
}


