fn main() {
    // When using system FFTW with threads, add the library search path.
    // On macOS with Homebrew, FFTW is typically in /opt/homebrew.
    if cfg!(feature = "fftw") {
        if let Ok(path) = std::process::Command::new("pkg-config")
            .args(["--libs-only-L", "fftw3"])
            .output()
        {
            let stdout = String::from_utf8_lossy(&path.stdout);
            for flag in stdout.split_whitespace() {
                if let Some(dir) = flag.strip_prefix("-L") {
                    println!("cargo:rustc-link-search=native={dir}");
                }
            }
        }
    }
}
