extern crate assert_cmd;
extern crate predicates;

use std::process;

use assert_cmd::prelude::*;
use predicates::prelude::*;

#[test]
fn resample() {
    let mut cmd = process::Command::cargo_example("testresample").unwrap();
    let expected = String::from("Quality: 10\n").into_bytes();
    cmd.assert().stdout(&predicate::eq(expected));
}
