#[cfg(doc)]
use crate::qag::Qag;

use crate::constants::{HeapItem, Myf64};
use ndarray::{array, Array1};
use std::collections::{BinaryHeap, HashMap};
/// Result of [integrate](Qag::integrate).
///
/// It contains the result [Array1], the error and optionally a [MoreInfo].
#[derive(Debug, Clone)]
pub struct QagIntegrationResult {
    pub result: Array1<f64>,
    pub abserr: f64,
    pub more_info: Option<MoreInfo>,
}

impl QagIntegrationResult {
    pub fn new_more_info(
        result: Array1<f64>,
        abserr: f64,
        neval: i32,
        last: usize,
        hash: HashMap<(Myf64, Myf64), Array1<f64>>,
        heap: BinaryHeap<HeapItem>,
    ) -> Self {
        Self {
            result,
            abserr,
            more_info: Some(MoreInfo::new(neval, last, hash, heap)),
        }
    }

    pub fn new(result: Array1<f64>, abserr: f64) -> Self {
        Self {
            result,
            abserr,
            more_info: None,
        }
    }

    pub fn new_error() -> Self {
        Self {
            result: array![0.0],
            abserr: 0.0,
            more_info: None,
        }
    }
}
/// Optional additional information result of [integrate](Qag::integrate).
///
/// It contains the number of function evaluation 'neval', the number of interval subdivision
/// 'last', the [HashMap] with the integration result for every sub-interval 'hash' and the [BinaryHeap]
/// with the error for every sub-interval 'heap'.
#[derive(Debug, Clone)]
pub struct MoreInfo {
    pub neval: i32,
    pub last: usize,
    pub hash: HashMap<(Myf64, Myf64), Array1<f64>>,
    pub heap: BinaryHeap<HeapItem>,
}

impl MoreInfo {
    pub fn new(
        neval: i32,
        last: usize,
        hash: HashMap<(Myf64, Myf64), Array1<f64>>,
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
