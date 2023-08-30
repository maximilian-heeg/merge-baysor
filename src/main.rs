use clap::Parser;
use std::process;

use merge_baysor::*;

fn main() {
    let args = Args::parse();

    if let Err(err) = run(args) {
        println!("{}", err);
        process::exit(1);
    }
}
