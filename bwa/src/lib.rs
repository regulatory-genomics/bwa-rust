use std::{ffi::{CStr, CString}, path::Path};
use noodles::sam::{self, header::record::value::map::ReferenceSequence};
use noodles::sam::header::record::value::Map;
use noodles::fastq;
use bstr;
use bstr::ByteSlice;
use itertools::Itertools;

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

    pub fn create_sam_header(&self) -> sam::Header {
        let ref_seqs = self.contig_names.iter().zip(self.contig_lengths.iter()).map(|(name, len)|
            (bstr::BString::from(name.as_str()), Map::<ReferenceSequence>::new(std::num::NonZeroUsize::try_from(*len).unwrap()))
        ).collect();
        sam::Header::builder().set_reference_sequences(ref_seqs).build()
    }
}

impl Drop for FMIndex {
    fn drop(&mut self) {
        unsafe {
            bwa_sys::bwa_idx_destroy(self.fm_index as *mut bwa_sys::bwaidx_t);
        }
    }
}

/// BWA opts object. Currently only default opts are enabled
#[derive(Debug, Copy, Clone)]
pub struct AlignerOpts {
    opts: bwa_sys::mem_opt_t,
}

impl Default for AlignerOpts {
    fn default() -> Self {
        Self::new()
    }
}

impl AlignerOpts {
    /// Create a `Bwaopts` object with default BWA parameters
    pub fn new() -> Self {
        let ptr = unsafe { bwa_sys::mem_opt_init() };
        let opts = unsafe { *ptr };
        unsafe { libc::free(ptr as *mut libc::c_void) };
        Self { opts }
    }

    pub fn get_actual_chunk_size(&self) -> usize {
        (self.opts.chunk_size * self.opts.n_threads).try_into().unwrap()
    }

    pub fn set_n_threads(mut self, n_threads: usize) -> Self {
        self.opts.n_threads = n_threads as i32;
        self
    }

    /// Set alignment scores
    pub fn set_scores(
        mut self,
        matchp: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> Self {
        self.opts.a = matchp;
        self.opts.b = mismatch;
        self.opts.o_del = gap_open;
        self.opts.o_ins = gap_open;
        self.opts.e_del = gap_extend;
        self.opts.e_ins = gap_extend;

        unsafe {
            bwa_sys::bwa_fill_scmat(matchp, mismatch, self.opts.mat.as_mut_ptr());
        }
        self
    }

    /// Set clipping score penalties
    pub fn set_clip_scores(mut self, clip5: i32, clip3: i32) -> Self {
        self.opts.pen_clip5 = clip5;
        self.opts.pen_clip3 = clip3;
        self
    }

    /// Set unpaired read penalty
    pub fn set_unpaired(mut self, unpaired: i32) -> Self {
        self.opts.pen_unpaired = unpaired;
        self
    }

    /// Mark shorter splits as secondary
    pub fn set_no_multi(mut self) -> Self {
        self.opts.flag |= bwa_sys::MEM_F_NO_MULTI as i32;
        self
    }

    fn pe_mode(mut self) -> Self {
        self.opts.flag |= bwa_sys::MEM_F_PE as i32;
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
    opts: AlignerOpts,
    header: sam::Header,
    pe_stats: PairedEndStats,
}

impl BurrowsWheelerAligner {
    pub fn new(
        index: FMIndex,
        opts: AlignerOpts,
        pe_stats: PairedEndStats,
    ) -> Self {
        let header = index.create_sam_header();
        Self { index, opts, header, pe_stats }
    }

    pub fn get_sam_header(&self) -> sam::Header {
        self.header.clone()
    }

    /// Return the chunk size (in terms of number of bases) used by the aligner
    pub fn chunk_size(&self) -> usize {
        self.opts.get_actual_chunk_size()
    }

    pub fn align_reads_iter<'a, I>(&'a self, mut records: I) -> impl Iterator<Item = sam::Record> + '_
    where
        I: Iterator<Item = fastq::Record> + 'a,
    {
        let max_chunk_length = self.opts.get_actual_chunk_size();
        std::iter::from_fn(move || {
            let mut chunk: Vec<_> = Vec::new();
            let mut accumulated_length = 0;

            // Take sequences until the accumulated length exceeds the max_chunk_length
            for record in records.by_ref() {
                accumulated_length += record.sequence().len();
                chunk.push(record);
                if accumulated_length >= max_chunk_length {
                    break;
                }
            }

            if chunk.is_empty() {
                None
            } else {
                let sams = self.align_reads(chunk.as_mut());
                Some(sams.into_iter())
            }
        }).flatten()
    }

    pub fn align_read_pairs_iter<'a, I>(&'a self, mut records: I) -> impl Iterator<Item = (sam::Record, sam::Record)> + '_
    where
        I: Iterator<Item = (fastq::Record, fastq::Record)> + 'a
    {
        let max_chunk_length = self.opts.get_actual_chunk_size();

        std::iter::from_fn(move || {
            let mut chunk: Vec<_> = Vec::new();
            let mut accumulated_length = 0;

            // Take sequences until the accumulated length exceeds the max_chunk_length
            for (record1, record2) in records.by_ref() {
                accumulated_length += record1.sequence().len() + record2.sequence().len();
                chunk.push(record1);
                chunk.push(record2);
                if accumulated_length >= max_chunk_length {
                    break;
                }
            }

            if chunk.is_empty() {
                None
            } else {
                let sams = self.align_reads( chunk.as_mut());
                Some(sams.into_iter())
            }
        }).flatten().tuples()
    }

    pub fn align_read_pairs(
        &self,
        records: &mut [(fastq::Record, fastq::Record)],
    ) -> impl ExactSizeIterator<Item = (sam::Record, sam::Record)> {
        let mut seqs = records.iter_mut().enumerate().flat_map(|(i, (fq1, fq2))|
            [new_bseq1_t(i, fq1), new_bseq1_t(i+1, fq2)]
        ).collect::<Vec<_>>();
        let sam = unsafe {
            let index = *self.index.fm_index;
            bwa_sys::mem_process_seqs(
                &self.opts.clone().pe_mode().opts, index.bwt, index.bns, index.pac,
                0,
                seqs.len().try_into().unwrap(),
                seqs.as_mut_ptr(), self.pe_stats.inner.as_ptr(),
            );
            seqs.into_iter().map(|seq| CStr::from_ptr(seq.sam).to_str().unwrap().as_bytes().try_into().unwrap())
                .tuples()
        };
        sam
    }

    pub fn align_reads(
        &self,
        records: &mut [fastq::Record],
    ) -> impl ExactSizeIterator<Item = sam::Record> {
        let mut seqs: Vec<_> = records.iter_mut().enumerate().map(|(i, fq)| new_bseq1_t(i, fq)).collect();
        let sam = unsafe {
            let index = *self.index.fm_index;
            bwa_sys::mem_process_seqs(
                &self.opts.opts, index.bwt, index.bns, index.pac,
                0,
                seqs.len().try_into().unwrap(),
                seqs.as_mut_ptr(), self.pe_stats.inner.as_ptr(),
            );
            seqs.into_iter().map(|seq| CStr::from_ptr(seq.sam).to_str().unwrap().as_bytes().try_into().unwrap())
        };
        sam
    }
}

fn new_bseq1_t(id: usize, fq: &mut fastq::Record) -> bwa_sys::bseq1_t {
    bwa_sys::bseq1_t {
        name: CString::new(fq.name().as_bytes()).unwrap().into_raw(),
        id: id.try_into().unwrap(),
        l_seq: fq.sequence().len() as i32,
        seq: fq.sequence_mut().as_mut_ptr() as *mut i8,
        qual: fq.quality_scores_mut().as_mut_ptr() as *mut i8,
        comment: std::ptr::null_mut(),
        sam: std::ptr::null_mut(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{fs, io::BufReader};
    use noodles::sam;
    use noodles::fastq;
    use flate2;

    #[test]
    fn test_align() {
        let fasta = "tests/data/ecoli.fa.gz";
        let location = "tests/temp";
        let _ = fs::remove_dir_all(location);
        let idx = FMIndex::new(fasta, location).unwrap();
        let opt = AlignerOpts::default();
        let aligner = BurrowsWheelerAligner::new(idx, opt, PairedEndStats::default());
        let header = aligner.get_sam_header();

        let mut writer = sam::io::Writer::new(std::fs::File::create("out.sam").unwrap());
        let reader = flate2::read::MultiGzDecoder::new(std::fs::File::open("tests/test.fq.gz").unwrap());
        let mut reader = fastq::Reader::new(BufReader::new(reader));
        aligner.align_reads_iter(reader.records().map(|x| x.unwrap())).for_each(|sam| {
            writer.write_record(&header, &sam).unwrap();
        });
    }

    #[test]
    fn test_align_pe() {
        let fasta = "tests/data/ecoli.fa.gz";
        let location = "tests/temp";
        let _ = fs::remove_dir_all(location);
        let idx = FMIndex::new(fasta, location).unwrap();
        let opt = AlignerOpts::default();
        let aligner = BurrowsWheelerAligner::new(idx, opt, PairedEndStats::default());
        let header = aligner.get_sam_header();

        let mut writer = sam::io::Writer::new(std::fs::File::create("out.sam").unwrap());
        let reader = flate2::read::MultiGzDecoder::new(std::fs::File::open("tests/test.fq.gz").unwrap());
        let mut reader = fastq::Reader::new(BufReader::new(reader));
        aligner.align_read_pairs_iter(reader.records().map(|x| x.unwrap()).tuples()).for_each(|(sam1, sam2)| {
            writer.write_record(&header, &sam1).unwrap();
            writer.write_record(&header, &sam2).unwrap();
        });
    }
}