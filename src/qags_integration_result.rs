use std::collections::{BinaryHeap, HashMap};
use crate::constants::{HeapItem, Myf64};

#[derive(Debug,Clone)]
pub struct QagsIntegrationResult<const N:usize> {
    pub result : [f64; N],
    pub abserr : f64,
    pub more_info : Option<MoreInfo<N>>,
}

impl<const N:usize> QagsIntegrationResult<N>{
    pub fn new_more_info(result : [f64; N], abserr : f64, neval : i32, last : usize,
               hash : HashMap<(Myf64,Myf64),[f64;N]>, heap : BinaryHeap<HeapItem>) -> Self {
        Self{
            result, abserr, more_info : Some(MoreInfo::new(neval,last,hash,heap))
        }
    }

    pub fn new(result : [f64; N], abserr : f64) -> Self {
        Self{
            result, abserr, more_info : None
        }
    }

    pub fn new_error() -> Self {
        Self{
            result : [0.0; N], abserr : 0.0, more_info : None
        }
    }
}

#[derive(Debug,Clone)]
pub struct MoreInfo<const N:usize>{
    pub neval : i32,
    pub last : usize,
    pub hash : HashMap<(Myf64,Myf64),[f64;N]>,
    pub heap : BinaryHeap<HeapItem>,
}

impl<const N:usize> MoreInfo<N>{
    pub fn new(neval : i32, last : usize, hash : HashMap<(Myf64,Myf64),[f64;N]>,
               heap : BinaryHeap<HeapItem>) -> Self {
        Self{
            neval, last, hash, heap
        }
    }

}