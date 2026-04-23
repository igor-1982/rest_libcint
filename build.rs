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

    // read CINT_DEV from environment variable;
    // if CINT_DEV set to "1", link directory from CINT_DIR
    if std::env::var("CINT_DEV").unwrap_or_default() == "1" {
        let cint_dir = std::env::var("CINT_DIR").expect("CINT_DIR environment variable must be set when CINT_DEV=1");
        println!("cargo:rustc-link-search=native={}", cint_dir);
        println!("cargo:rustc-link-lib=cint");
    }
}
