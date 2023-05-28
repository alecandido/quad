mod constants;
mod qag_vec_norm_integrator_result;
mod qage_vec_norm2;
mod qage_vec_norm3;
mod qage_vec_norm_parall;
mod qk61_vec_norm2;
mod result_state;
mod qk;
mod qk61;
mod qk51;
mod qk41;
mod qk31;
mod qk21;
mod qk15;
mod qag_vec_norm_integration_result;
mod qag;
mod qagv2;
mod semi_infinite_function;
mod qags;
mod qags_integration_result;
mod qags_integrator_result;
mod qags2;

use pyo3::prelude::*;
use pyo3::types::PyTuple;
use crate::constants::{Myf64, norm_vec};
use crate::qag::Qag;
use crate::qag_vec_norm_integration_result::QagVecNormIntegrationResult;
use crate::qag_vec_norm_integrator_result::QagVecNormIntegratorResult;
use crate::qage_vec_norm3::QagVecNorm3;
use crate::qags2::Qags2;
use crate::qags::Qags;
use crate::qags_integration_result::MoreInfo;
use crate::qagv2::Qagv2;
use crate::result_state::ResultState::MaxIteration;


#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

fn lambda_eval<const n : usize>(ob : &PyAny, z : f64) -> [f64;n]{
    let f = |x:f64| {
        let y = (x,);
        ob.call1(y).unwrap().extract::<[f64;n]>().unwrap()
    };
    f(z)
}


#[pyfunction]
fn qag(ob : &PyAny, a: f64, b: f64, epsabs : f64, epsrel : f64, key : i32, dim : i32)
    -> Py<PyTuple> {

    let res = Python::with_gil(|py| {
        let qag = Qag { key, limit: 1000000 };

        match dim {
            1 => {
                let f = |x: f64| lambda_eval::<1>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },
            2 => {
                let f = |x: f64| lambda_eval::<2>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },

            3 => {
                let f = |x: f64| lambda_eval::<3>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },

            4 => {
                let f = |x: f64| lambda_eval::<4>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },

            5 => {
                let f = |x: f64| lambda_eval::<5>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },

            6 => {
                let f = |x: f64| lambda_eval::<6>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },

            7 => {
                let f = |x: f64| lambda_eval::<7>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },

            8 => {
                let f = |x: f64| lambda_eval::<8>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },


            _ => {
                let f = |x: f64| lambda_eval::<4>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },
        }
    });
    res





}

#[pyfunction]
fn qagv2(ob : &PyAny, a: f64, b: f64, epsabs : f64, epsrel : f64, key : i32, dim : i32)
       -> Py<PyTuple> {

    let res = Python::with_gil(|py| {
        let qag = Qagv2 { key, limit: 1000000 };

        match dim {
            1 => {
                let f = |x: f64| lambda_eval::<1>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },
            2 => {
                let f = |x: f64| lambda_eval::<2>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },

            3 => {
                let f = |x: f64| lambda_eval::<3>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },

            4 => {
                let f = |x: f64| lambda_eval::<4>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },

            5 => {
                let f = |x: f64| lambda_eval::<5>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },

            6 => {
                let f = |x: f64| lambda_eval::<6>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },

            7 => {
                let f = |x: f64| lambda_eval::<7>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },

            8 => {
                let f = |x: f64| lambda_eval::<8>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },


            _ => {
                let f = |x: f64| lambda_eval::<4>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel);
                let mut iter = res.integration_result.result.to_vec();
                iter.push(res.integration_result.abserr);
                let a : Py<PyTuple> = PyTuple::new(py,iter.into_iter()).into();
                a
            },
        }
    });
    res
}


#[pyfunction]
fn qagv3(ob : &PyAny, a: f64, b: f64, epsabs : f64, epsrel : f64, key : i32, dim : i32)
         -> Py<MyClass> {

    let out = Python::with_gil(|py| {
        let qag = Qagv2 { key, limit: 1000000 };

        match dim {
            1 => {
                let f = |x: f64| lambda_eval::<1>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },
            2 => {
                let f = |x: f64| lambda_eval::<2>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },

            3 => {
                let f = |x: f64| lambda_eval::<3>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },

            4 => {
                let f = |x: f64| lambda_eval::<4>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },

            5 => {
                let f = |x: f64| lambda_eval::<5>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },

            6 => {
                let f = |x: f64| lambda_eval::<6>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },

            7 => {
                let f = |x: f64| lambda_eval::<7>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },

            8 => {
                let f = |x: f64| lambda_eval::<8>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },


            _ => {
                let f = |x: f64| lambda_eval::<4>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },
        }
    });
    out
}

#[pyfunction]
fn qags(ob : &PyAny, a: f64, b: f64, epsabs : Option<f64>, epsrel : Option<f64>, key : Option<i32>, limit : Option<usize>, dim : i32,
        points : Option<Vec<f64>>)
         -> Py<MyClass> {

    let pointss = {
    if points.is_some() {
        points.unwrap()
    }
    else {
        [0.0;0].to_vec()
    }
    };
    let limitt = {
        if limit.is_some(){
            limit.unwrap()
        }
        else { 1000 }
    };
    let keyy = {
        if key.is_some(){
            key.unwrap()
        }
        else{
            2
        }
    };
    let epsabss = {
        if epsabs.is_some(){
            epsabs.unwrap()
        }
        else {
            1e-4
        }
    };
    let epsrell = {
        if epsrel.is_some(){
            epsrel.unwrap()
        }
        else {
            1e-4
        }
    };



    let out = Python::with_gil(|py| {
        let qag = Qags{key : keyy,limit : limitt,points : pointss};

        match dim {
            1 => {
                let f = |x: f64| lambda_eval::<1>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },
            2 => {
                let f = |x: f64| lambda_eval::<2>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },

            3 => {
                let f = |x: f64| lambda_eval::<3>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },

            4 => {
                let f = |x: f64| lambda_eval::<4>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },

            5 => {
                let f = |x: f64| lambda_eval::<5>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },

            6 => {
                let f = |x: f64| lambda_eval::<6>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },

            7 => {
                let f = |x: f64| lambda_eval::<7>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },

            8 => {
                let f = |x: f64| lambda_eval::<8>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },


            _ => {
                let f = |x: f64| lambda_eval::<4>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,neval,last) = (res.result.to_vec(),res.abserr,res.neval,res.last);
                Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
            },
        }
    });
    out
}


#[pyfunction]
fn qags2(ob : &PyAny, a: f64, b: f64, epsabs : Option<f64>, epsrel : Option<f64>, key : Option<i32>, limit : Option<usize>, dim : i32,
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
        if limit.is_some(){
            limit.unwrap()
        }
        else { 1000 }
    };
    let keyy = {
        if key.is_some(){
            key.unwrap()
        }
        else{
            2
        }
    };
    let epsabss = {
        if epsabs.is_some(){
            epsabs.unwrap()
        }
        else {
            1e-4
        }
    };
    let epsrell = {
        if epsrel.is_some(){
            epsrel.unwrap()
        }
        else {
            1e-4
        }
    };
    let more_infoo = {
        if more_info.is_some() {
            more_info.unwrap()
        }
        else { false }
    };



    let out = Python::with_gil(|py| {
        let qag = Qags2{key : keyy,limit : limitt, points : pointss, more_info : more_infoo};

        match dim {
            1 => {
                let f = |x: f64| lambda_eval::<1>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,more_inf) = (res.result.to_vec(),res.abserr,res.more_info);
                if more_inf.is_none(){
                    Py::new(py, QagsResult {result,abserr,more_info : None}).unwrap()
                }
                else {
                    let mut more_inf_py : Vec<(f64,f64,f64,Vec<f64>)> = vec![];
                    let mut more_inf_unwrapped = more_inf.unwrap();
                    let (mut hash, mut heap) = (more_inf_unwrapped.hash, more_inf_unwrapped.heap);
                    let (neval, last) = (more_inf_unwrapped.neval, more_inf_unwrapped.last);
                    for _k in 0..heap.len(){
                        let old_interval = heap.pop().unwrap();
                        let ((x,y),old_err) = (old_interval.interval,old_interval.err);
                        let old_res = hash.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
                        more_inf_py.push((x,y,old_err,old_res.to_vec()));
                    }
                    Py::new(py, QagsResult {result,abserr,more_info : Some((neval,last,more_inf_py))}).unwrap()
                }

            },
            2 => {
                let f = |x: f64| lambda_eval::<2>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,more_inf) = (res.result.to_vec(),res.abserr,res.more_info);
                if more_inf.is_none(){
                    Py::new(py, QagsResult {result,abserr,more_info : None}).unwrap()
                }
                else {
                    let mut more_inf_py : Vec<(f64,f64,f64,Vec<f64>)> = vec![];
                    let mut more_inf_unwrapped = more_inf.unwrap();
                    let (mut hash, mut heap) = (more_inf_unwrapped.hash, more_inf_unwrapped.heap);
                    let (neval, last) = (more_inf_unwrapped.neval, more_inf_unwrapped.last);
                    for _k in 0..heap.len(){
                        let old_interval = heap.pop().unwrap();
                        let ((x,y),old_err) = (old_interval.interval,old_interval.err);
                        let old_res = hash.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
                        more_inf_py.push((x,y,old_err,old_res.to_vec()));
                    }
                    Py::new(py, QagsResult {result,abserr,more_info : Some((neval,last,more_inf_py))}).unwrap()
                }
            },

            3 => {
                let f = |x: f64| lambda_eval::<3>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,more_inf) = (res.result.to_vec(),res.abserr,res.more_info);
                if more_inf.is_none(){
                    Py::new(py, QagsResult {result,abserr,more_info : None}).unwrap()
                }
                else {
                    let mut more_inf_py : Vec<(f64,f64,f64,Vec<f64>)> = vec![];
                    let mut more_inf_unwrapped = more_inf.unwrap();
                    let (mut hash, mut heap) = (more_inf_unwrapped.hash, more_inf_unwrapped.heap);
                    let (neval, last) = (more_inf_unwrapped.neval, more_inf_unwrapped.last);
                    for _k in 0..heap.len(){
                        let old_interval = heap.pop().unwrap();
                        let ((x,y),old_err) = (old_interval.interval,old_interval.err);
                        let old_res = hash.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
                        more_inf_py.push((x,y,old_err,old_res.to_vec()));
                    }
                    Py::new(py, QagsResult {result,abserr,more_info : Some((neval,last,more_inf_py))}).unwrap()
                }
            },

            4 => {
                let f = |x: f64| lambda_eval::<4>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,more_inf) = (res.result.to_vec(),res.abserr,res.more_info);
                if more_inf.is_none(){
                    Py::new(py, QagsResult {result,abserr,more_info : None}).unwrap()
                }
                else {
                    let mut more_inf_py : Vec<(f64,f64,f64,Vec<f64>)> = vec![];
                    let mut more_inf_unwrapped = more_inf.unwrap();
                    let (mut hash, mut heap) = (more_inf_unwrapped.hash, more_inf_unwrapped.heap);
                    let (neval, last) = (more_inf_unwrapped.neval, more_inf_unwrapped.last);
                    for _k in 0..heap.len(){
                        let old_interval = heap.pop().unwrap();
                        let ((x,y),old_err) = (old_interval.interval,old_interval.err);
                        let old_res = hash.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
                        more_inf_py.push((x,y,old_err,old_res.to_vec()));
                    }
                    Py::new(py, QagsResult {result,abserr,more_info : Some((neval,last,more_inf_py))}).unwrap()
                }
            },

            5 => {
                let f = |x: f64| lambda_eval::<5>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,more_inf) = (res.result.to_vec(),res.abserr,res.more_info);
                if more_inf.is_none(){
                    Py::new(py, QagsResult {result,abserr,more_info : None}).unwrap()
                }
                else {
                    let mut more_inf_py : Vec<(f64,f64,f64,Vec<f64>)> = vec![];
                    let mut more_inf_unwrapped = more_inf.unwrap();
                    let (mut hash, mut heap) = (more_inf_unwrapped.hash, more_inf_unwrapped.heap);
                    let (neval, last) = (more_inf_unwrapped.neval, more_inf_unwrapped.last);
                    for _k in 0..heap.len(){
                        let old_interval = heap.pop().unwrap();
                        let ((x,y),old_err) = (old_interval.interval,old_interval.err);
                        let old_res = hash.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
                        more_inf_py.push((x,y,old_err,old_res.to_vec()));
                    }
                    Py::new(py, QagsResult {result,abserr,more_info : Some((neval,last,more_inf_py))}).unwrap()
                }
            },

            6 => {
                let f = |x: f64| lambda_eval::<6>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,more_inf) = (res.result.to_vec(),res.abserr,res.more_info);
                if more_inf.is_none(){
                    Py::new(py, QagsResult {result,abserr,more_info : None}).unwrap()
                }
                else {
                    let mut more_inf_py : Vec<(f64,f64,f64,Vec<f64>)> = vec![];
                    let mut more_inf_unwrapped = more_inf.unwrap();
                    let (mut hash, mut heap) = (more_inf_unwrapped.hash, more_inf_unwrapped.heap);
                    let (neval, last) = (more_inf_unwrapped.neval, more_inf_unwrapped.last);
                    for _k in 0..heap.len(){
                        let old_interval = heap.pop().unwrap();
                        let ((x,y),old_err) = (old_interval.interval,old_interval.err);
                        let old_res = hash.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
                        more_inf_py.push((x,y,old_err,old_res.to_vec()));
                    }
                    Py::new(py, QagsResult {result,abserr,more_info : Some((neval,last,more_inf_py))}).unwrap()
                }
            },

            7 => {
                let f = |x: f64| lambda_eval::<7>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,more_inf) = (res.result.to_vec(),res.abserr,res.more_info);
                if more_inf.is_none(){
                    Py::new(py, QagsResult {result,abserr,more_info : None}).unwrap()
                }
                else {
                    let mut more_inf_py : Vec<(f64,f64,f64,Vec<f64>)> = vec![];
                    let mut more_inf_unwrapped = more_inf.unwrap();
                    let (mut hash, mut heap) = (more_inf_unwrapped.hash, more_inf_unwrapped.heap);
                    let (neval, last) = (more_inf_unwrapped.neval, more_inf_unwrapped.last);
                    for _k in 0..heap.len(){
                        let old_interval = heap.pop().unwrap();
                        let ((x,y),old_err) = (old_interval.interval,old_interval.err);
                        let old_res = hash.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
                        more_inf_py.push((x,y,old_err,old_res.to_vec()));
                    }
                    Py::new(py, QagsResult {result,abserr,more_info : Some((neval,last,more_inf_py))}).unwrap()
                }
            },

            8 => {
                let f = |x: f64| lambda_eval::<8>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,more_inf) = (res.result.to_vec(),res.abserr,res.more_info);
                if more_inf.is_none(){
                    Py::new(py, QagsResult {result,abserr,more_info : None}).unwrap()
                }
                else {
                    let mut more_inf_py : Vec<(f64,f64,f64,Vec<f64>)> = vec![];
                    let mut more_inf_unwrapped = more_inf.unwrap();
                    let (mut hash, mut heap) = (more_inf_unwrapped.hash, more_inf_unwrapped.heap);
                    let (neval, last) = (more_inf_unwrapped.neval, more_inf_unwrapped.last);
                    for _k in 0..heap.len(){
                        let old_interval = heap.pop().unwrap();
                        let ((x,y),old_err) = (old_interval.interval,old_interval.err);
                        let old_res = hash.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
                        more_inf_py.push((x,y,old_err,old_res.to_vec()));
                    }
                    Py::new(py, QagsResult {result,abserr,more_info : Some((neval,last,more_inf_py))}).unwrap()
                }
            },

            _ => {
                let f = |x: f64| lambda_eval::<1>(ob, x);
                let res = qag.qintegrate(&f, a, b, epsabss, epsrell).unwrap();
                let (result,abserr,more_inf) = (res.result.to_vec(),res.abserr,res.more_info);
                if more_inf.is_none(){
                    Py::new(py, QagsResult {result,abserr,more_info : None}).unwrap()
                }
                else {
                    let mut more_inf_py : Vec<(f64,f64,f64,Vec<f64>)> = vec![];
                    let mut more_inf_unwrapped = more_inf.unwrap();
                    let (mut hash, mut heap) = (more_inf_unwrapped.hash, more_inf_unwrapped.heap);
                    let (neval, last) = (more_inf_unwrapped.neval, more_inf_unwrapped.last);
                    for _k in 0..heap.len(){
                        let old_interval = heap.pop().unwrap();
                        let ((x,y),old_err) = (old_interval.interval,old_interval.err);
                        let old_res = hash.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
                        more_inf_py.push((x,y,old_err,old_res.to_vec()));
                    }
                    Py::new(py, QagsResult {result,abserr,more_info : Some((neval,last,more_inf_py))}).unwrap()
                }

            },
        }
    });
    out
}


#[pyfunction]
fn test_my_class() -> Py<MyClass>{
    let res = Python::with_gil(|py| {
        let result = vec![0.0,2.0];
        let abserr = 0.2;
        let neval = 1;
        let last = 10;

        Py::new(py,MyClass{result,abserr,neval,last}).unwrap()
    });
    res
}

#[pyclass]
struct MyClass {
    #[pyo3(get, set)]
    pub result : Vec<f64>,
    #[pyo3(get, set)]
    pub abserr : f64,
    #[pyo3(get, set)]
    pub neval : i32,
    #[pyo3(get, set)]
    pub last : usize,
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
    pub neval : i32,
    pub last : usize,
    pub hash : Vec<(f64,f64,f64,Vec<f64>)>,
}


#[pyfunction]
fn test_option(list : Option<Vec<f64>>) -> () {
    if list.is_some() {
        let vec = list.unwrap().clone();
        let qag = Qags { key: 6, limit: 1000, points: vec };
        println!("{:?}", qag);
    } else {
        let vec = [0.0; 0].to_vec();
        let qag = Qags { key: 6, limit: 1000, points: vec };
        println!("{:?}", qag);
    }
}


/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn quad(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(qag, m)?)?;
    m.add_function(wrap_pyfunction!(qagv2, m)?)?;
    m.add_function(wrap_pyfunction!(qagv3, m)?)?;
    m.add_function(wrap_pyfunction!(qags, m)?)?;
    m.add_function(wrap_pyfunction!(qags2, m)?)?;
    m.add_function(wrap_pyfunction!(test_my_class, m)?)?;
    m.add_function(wrap_pyfunction!(test_option, m)?)?;
    m.add_class::<MyClass>()?;
    Ok(())
}






#[cfg(test)]
mod tests {
    use std::time::Instant;
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


    #[test]
    fn array_vs_vec(){

        let start = Instant::now();
        let mut array = [0.0;128];
        for k in 0..128{
            array[k] = k as f64;
        }
        println!("array time : {:?}",start.elapsed());

        let start = Instant::now();
        let mut vec = vec![];
        for k in 0..15{
            vec.push(k as f64);
        }
        println!("vec time : {:?}",start.elapsed());


    }




}


