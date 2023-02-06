use crate::*;
use std::sync::{Arc, Mutex};

pub struct FunctVector {
    pub components: Vec<Box<dyn Function + Send + Sync >>,
}

pub struct FunVector {
    pub components : Vec<Arc<Mutex<dyn Function + Send + Sync>>>,
}

