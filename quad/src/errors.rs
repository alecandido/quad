use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub enum QagError {
    Invalid,
    MaxIteration,
    BadTolerance,
    BadFunction,
    Diverge,
}

impl fmt::Display for QagError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let error_message : &str;
        match self{
            QagError::Invalid=> error_message = INVALID_ERROR_MESSAGE,
            QagError::MaxIteration => error_message = MAX_ITERATION_ERROR_MESSAGE,
            QagError::BadTolerance=> error_message = BAD_TOLERANCE_ERROR_MESSAGE,
            QagError::BadFunction=> error_message = BAD_FUNCTION_ERROR_MESSAGE,
            QagError::Diverge=> error_message = DIVERGE_ERROR_MESSAGE,
        }
        write!(f, "{}", error_message)
    }
}



pub const MAX_ITERATION_ERROR_MESSAGE: &str =
    "Maximum number of subdivisions allowed has been achieved.
            One can allow more subdivisions by increasing the value of limit. However, if this
            yields no improvement it is rather advised to analyze the integrand in order to
            determine the integration difficulties. If the position of a local difficulty can be
            determined(e.g. singularity, discontinuity within the interval) one will probably gain
            from splitting up the interval at this point and calling the integrator on the
            subranges. If possible, an appropriate special-purpose integrator should be used
            which is designed for handling the type of difficulty involved.";
pub const BAD_TOLERANCE_ERROR_MESSAGE: &str =
    "The occurrence of roundoff error is detected, which \
            prevents the requested tolerance from being achieved.";
pub const INVALID_ERROR_MESSAGE: &str = "The input is invalid, because epsabs <= 0 and \
            epsrel < max(50 * rel.mach.acc.,0.5d-28)";
pub const BAD_FUNCTION_ERROR_MESSAGE: &str = "Extremely bad integrand behaviour occurs at some \
            points of the integration interval.";
pub const DIVERGE_ERROR_MESSAGE: &str = "The integral is probably divergent, or slowly convergent.\
            It must be noted that divergence can occur with any other value of ResultState.";
