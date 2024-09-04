use std::{ffi::{CStr, CString}, path::Path};
use noodles::sam;

#[derive(Debug)]
pub struct IndexError(String);

/// A BWA reference object to perform alignments to.
/// Must be loaded from a BWA index created with `bwa index`
pub struct FMIndex {
    fm_index: *const bwa_sys::bwaidx_t,
    contig_names: Vec<String>,
    contig_lengths: Vec<usize>,
}

impl FMIndex {
    pub fn new<P1: AsRef<Path>, P2: AsRef<Path>>(fasta: P1, location: P2) -> Result<Self, IndexError> {
        let fasta_file = CString::new(fasta.as_ref().to_str().unwrap()).unwrap();
        let prefix = CString::new(location.as_ref().to_str().unwrap()).unwrap();
        unsafe {
            bwa_sys::bwa_idx_build(
                fasta_file.as_ptr(), 
                prefix.as_ptr(), 
                bwa_sys::BWTALGO_BWTSW.try_into().unwrap(),
                10000000,
            );
        }
        Self::read(location)
    }

    /// Load a BWA reference from disk. Pass the fasta filename of the
    /// original reference as `path`
    pub fn read<P: AsRef<Path>>(path: P) -> Result<Self, IndexError> {
        let idx_file = CString::new(path.as_ref().to_str().unwrap()).unwrap();
        let idx = unsafe {
            bwa_sys::bwa_idx_load(idx_file.as_ptr(), bwa_sys::BWA_IDX_ALL.try_into().unwrap())
        };

        let mut contig_names = Vec::new();
        let mut contig_lengths = Vec::new();
        let num_contigs = unsafe { (*(*idx).bns).n_seqs };

        for i in 0..num_contigs as isize {
            unsafe {
                let name = CStr::from_ptr((*(*(*idx).bns).anns.offset(i)).name);
                let sz = (*(*(*idx).bns).anns.offset(i)).len;

                let name_string = name.to_owned().into_string().unwrap();
                contig_names.push(name_string);
                contig_lengths.push(sz as usize)
            }
        }

        Ok(Self { fm_index: idx, contig_names, contig_lengths })
    }
}

impl Drop for FMIndex {
    fn drop(&mut self) {
        unsafe {
            bwa_sys::bwa_idx_destroy(self.fm_index as *mut bwa_sys::bwaidx_t);
        }
    }
}

/// BWA settings object. Currently only default settings are enabled
pub struct AlignerSettings {
    settings: bwa_sys::mem_opt_t,
}

impl Default for AlignerSettings {
    fn default() -> Self {
        Self::new()
    }
}

impl AlignerSettings {
    /// Create a `BwaSettings` object with default BWA parameters
    pub fn new() -> Self {
        let ptr = unsafe { bwa_sys::mem_opt_init() };
        let settings = unsafe { *ptr };
        unsafe { libc::free(ptr as *mut libc::c_void) };
        Self { settings }
    }

    /// Set alignment scores
    pub fn set_scores(
        mut self,
        matchp: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> Self {
        self.settings.a = matchp;
        self.settings.b = mismatch;
        self.settings.o_del = gap_open;
        self.settings.o_ins = gap_open;
        self.settings.e_del = gap_extend;
        self.settings.e_ins = gap_extend;

        unsafe {
            bwa_sys::bwa_fill_scmat(matchp, mismatch, self.settings.mat.as_mut_ptr());
        }
        self
    }

    /// Set clipping score penalties
    pub fn set_clip_scores(mut self, clip5: i32, clip3: i32) -> Self {
        self.settings.pen_clip5 = clip5;
        self.settings.pen_clip3 = clip3;
        self
    }

    /// Set unpaired read penalty
    pub fn set_unpaired(mut self, unpaired: i32) -> Self {
        self.settings.pen_unpaired = unpaired;
        self
    }

    /// Mark shorter splits as secondary
    pub fn set_no_multi(mut self) -> Self {
        self.settings.flag |= bwa_sys::MEM_F_NO_MULTI as i32;
        self
    }
}

/// Paired-end statistics structure used by BWA to score paired-end reads
pub struct PairedEndStats {
    inner: [bwa_sys::mem_pestat_t; 4],
}

impl Default for PairedEndStats {
    fn default() -> Self {
        Self::simple(200.0, 100.0, 35, 600)
    }
}

impl PairedEndStats {
    /// Generate a 'simple' paired-end read structure that standard forward-reverse
    /// pairs as created by TruSeq, Nextera, or Chromium Genome sample preparations.
    pub fn simple(avg: f64, std: f64, low: i32, high: i32) -> PairedEndStats {
        let pe_stat_null = || bwa_sys::mem_pestat_t {
            failed: 1,
            low: 0,
            high: 0,
            avg: 0.0,
            std: 100.0,
        };

        let pes = [
            pe_stat_null(),
            bwa_sys::mem_pestat_t {
                failed: 0,
                low,
                high,
                avg,
                std,
            },
            pe_stat_null(),
            pe_stat_null(),
        ];

        PairedEndStats { inner: pes }
    }
}

/// A BWA aligner. Carries everything required to align
/// reads to a reference and generate BAM records.
pub struct BurrowsWheelerAligner {
    index: FMIndex,
    settings: AlignerSettings,
    pe_stats: PairedEndStats,
}

impl BurrowsWheelerAligner {
    pub fn new(
        index: FMIndex,
        settings: AlignerSettings,
        pe_stats: PairedEndStats,
    ) -> Self {
        Self {
            index,
            settings,
            pe_stats,
        }
    }

    /// Align a read-pair to the reference.
    pub fn align_read(&self, name: &[u8], seq: &[u8], qual: &[u8]) -> sam::Record
    {
        unsafe {
            let sam_ptr = bwa_sys::mem_align(self.index.fm_index, name.as_ptr() as *mut i8, seq.as_ptr() as *mut i8, qual.as_ptr() as *mut i8);
            let sam = CStr::from_ptr(sam_ptr).to_str().unwrap().as_bytes().try_into().unwrap();
            libc::free(sam_ptr as *mut libc::c_void);
            sam
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use noodles::sam;

    #[test]
    fn test_index() {
        let data = [
            ( b"@chr_1561275_1561756_1:0:0_2:0:0_5c/1".to_vec(),
              b"GCATCGATAAGCAGGTCAAATTCTCCCGTCATTATCACCTCTGCTACTTAAATTTCCCGCTTTATAAGCCGATTACGGCCTGGCATTACCCTATCCATAATTTAGGTGGGATGCCCGGTGCGTGGTTGGCAGATCCGCTGTTCTTTATTT".to_vec(),
              b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222".to_vec(),
            ),
            ( b"@chr_727436_727956_3:0:0_1:0:0_0/1".to_vec(),
              b"GATGGCTGCGCAAGGGTTCTTACTGATCGCCACGTTTTTACTGGTGTTAATGGTGCTGGCGCGTCCTTTAGGCAGCGGGCTGGCGCGGCTGATTAATGACATTCCTCTTCCCGGTACAACGGGCGTTGAGCGCGAACTTTTTCGCGCACT".to_vec(),
              b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222".to_vec(),
            ),
        ];

        let fasta = "tests/data/ecoli.fa.gz";
        let location = "tests/temp";
        let _ = fs::remove_dir_all(location);
        let idx = FMIndex::new(fasta, location).unwrap();
        let aligner = BurrowsWheelerAligner::new(idx, AlignerSettings::default(), PairedEndStats::default());

        //let mut writer = sam::io::Writer::new(Vec::new());
        let sam = aligner.align_read(&data[0].0, &data[0].1, &data[0].2);
        println!("{:?}", sam);
    }
}