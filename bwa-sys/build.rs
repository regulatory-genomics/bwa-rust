// make -C bwa-sys/bwa/ -n libbwa.a | grep -o -E "[A-Za-z0-9_]+\.c"
const FILES: &[&str] = &[
    "ext/bwa/utils.c",
    "ext/bwa/QSufSort.c",
    "ext/bwa/kthread.c",
    "ext/bwa/kstring.c",
    "ext/bwa/rope.c",
    "ext/bwa/is.c",
    "ext/bwa/rle.c",
    "ext/bwa/ksw.c",
    "ext/bwa/bwt.c",
    "ext/bwa/bntseq.c",
    "ext/bwa/bwa.c",
    "ext/bwa/bwamem.c",
    "ext/bwa/bwtindex.c",
    "ext/bwa/bwamem_pair.c",
    "ext/bwa/bwamem_extra.c",
    "ext/bwa/malloc_wrap.c",
    "ext/bwa/bwt_gen.c",
];

fn main() {
    for file in FILES {
        println!("cargo:rerun-if-changed={}", file);
    }
    cc::Build::new()
        .warnings(false)
        .extra_warnings(false)
        .flag("-O3")
        .flag("-fPIC")
        .files(FILES)
        .compile("libbwa.a");

    println!("cargo:rustc-link-lib=z");
}