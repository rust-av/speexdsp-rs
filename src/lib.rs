extern crate speexdsp_sys;

pub use speexdsp_sys::*;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
