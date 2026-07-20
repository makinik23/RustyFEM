use clap::Parser;

#[derive(Debug, Parser)]
#[command(name = "rusty-fem", version, about = "Educational finite element method solver written in Rust")]
struct Cli;

fn main() {
    Cli::parse();
    println!("RustyFEM project scaffold is ready. The 1D bar solver is not implemented yet.");
}
