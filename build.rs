fn build_ecp() {
    cc::Build::new()
        .file("src/cecp/f2c_dgemm.c")
        .file("src/cecp/nr_ecp.c")
        .file("src/cecp/nr_ecp_deriv.c")
        .flag_if_supported("-Wno-unused-parameter")
        .flag_if_supported("-Wno-implicit-function-declaration")
        .flag_if_supported("-Wno-parentheses")
        .compile("cecp");
    println!("cargo::rerun-if-changed=src/cecp");
}

fn main() {
    build_ecp();
}
