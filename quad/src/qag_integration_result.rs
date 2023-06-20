use crate::constants::{HeapItem, Myf64};
use std::collections::{BinaryHeap, HashMap};

#[derive(Debug, Clone)]
pub struct QagIntegrationResult {
    pub result: Vec<f64>,
    pub abserr: f64,
    pub more_info: Option<MoreInfo>,
}

impl QagIntegrationResult {
    pub fn new_more_info(
        result: Vec<f64>,
        abserr: f64,
        neval: i32,
        last: usize,
        hash: HashMap<(Myf64, Myf64), Vec<f64>>,
        heap: BinaryHeap<HeapItem>,
    ) -> Self {
        Self {
            result,
            abserr,
            more_info: Some(MoreInfo::new(neval, last, hash, heap)),
        }
    }

    pub fn new(result: Vec<f64>, abserr: f64) -> Self {
        Self {
            result,
            abserr,
            more_info: None,
        }
    }

    pub fn new_error() -> Self {
        Self {
            result: vec![0.0],
            abserr: 0.0,
            more_info: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct MoreInfo {
    pub neval: i32,
    pub last: usize,
    pub hash: HashMap<(Myf64, Myf64), Vec<f64>>,
    pub heap: BinaryHeap<HeapItem>,
}

impl MoreInfo {
    pub fn new(
        neval: i32,
        last: usize,
        hash: HashMap<(Myf64, Myf64), Vec<f64>>,
        heap: BinaryHeap<HeapItem>,
    ) -> Self {
        Self {
            neval,
            last,
            hash,
            heap,
        }
    }
}
