use crate::*;
use std::sync::{Arc, Mutex};


#[derive(Clone)]
pub struct FnPa {
    pub components : Arc<dyn Fn(f64)->f64 + Send + Sync>,
}

