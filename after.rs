use quickdna::{NucleotideIter, NucleotideLike, ToNucleotideLike, TranslationTable};

/// Provides iterator of all windows (both nucleotides and amino acids) of supplied DNA.
#[derive(Clone, Debug, Default)]
pub struct DnaAndProteinWindows {
    dna_windows: DnaWindows,
    protein_windows: ProteinWindows,
}

impl DnaAndProteinWindows {
    pub fn from_dna<D, I, N>(dna: D) -> Self
    where
        D: IntoIterator<IntoIter = I>,
        I: DoubleEndedIterator<Item = N> + ExactSizeIterator + Clone,
        N: ToNucleotideLike,
    {
        let dna = dna.into_iter();
        Self {
            dna_windows: DnaWindows::from_dna(dna.clone()),
            protein_windows: ProteinWindows::from_dna(dna),
        }
    }

    /// Returns an iterator over the windows and their original position in the full FASTA
    ///
    /// The indexes returned by `enumerate_windows` are not strictly
    /// monotonically increasing because the current implementation (e.g.)
    /// yields all forward DNA windows before all RC DNA windows.
    pub fn enumerate_windows(
        &self,
        dna_window_len: usize,
        aa_window_len: usize,
    ) -> impl ExactSizeIterator<Item = (usize, &str)> {
        let dna_windows = self.dna_windows.enumerate_windows(dna_window_len);
        let protein_windows = self.protein_windows.enumerate_windows(aa_window_len);
        WithExactSize(dna_windows.chain(protein_windows))
    }
}

/// Provides iterator of nucleotide windows of supplied DNA.
#[derive(Clone, Debug, Default)]
pub struct DnaWindows {
    dna: String,
    dna_rc: String,
}

impl DnaWindows {
    pub fn from_dna<D, I, N>(dna: D) -> Self
    where
        D: IntoIterator<IntoIter = I>,
        I: DoubleEndedIterator<Item = N> + Clone,
        N: ToNucleotideLike,
    {
        let dna = dna.into_iter();

        let dna_rc = dna
            .clone()
            .reverse_complement()
            .map(|n| n.to_nucleotide_like().to_ascii())
            .collect();
        let dna_rc = String::from_utf8(dna_rc).unwrap();

        let dna = dna.map(|n| n.to_nucleotide_like().to_ascii()).collect();
        let dna = String::from_utf8(dna).unwrap();

        Self { dna, dna_rc }
    }

    /// Returns an iterator over the windows and their original position in the full FASTA
    ///
    /// The indexes returned by `enumerate_windows` are not strictly
    /// monotonically increasing because the current implementation (e.g.)
    /// yields all forward windows before all RC windows.
    pub fn enumerate_windows(
        &self,
        window_len: usize,
    ) -> impl ExactSizeIterator<Item = (usize, &str)> {
        let forward_windows = ascii_str_windows(&self.dna, window_len).enumerate();
        let rc_windows = ascii_str_windows(&self.dna_rc, window_len)
            .rev()
            .enumerate();
        WithExactSize(forward_windows.chain(rc_windows))
    }
}

/// Provides iterator of amino acids windows of supplied DNA.
#[derive(Clone, Debug, Default)]
pub struct ProteinWindows {
    aas: [String; 3],
    aas_rc: [String; 3],
}

impl ProteinWindows {
    pub fn from_dna<D, I, N>(dna: D) -> Self
    where
        D: IntoIterator<IntoIter = I>,
        I: DoubleEndedIterator<Item = N> + ExactSizeIterator + Clone,
        N: ToNucleotideLike,
    {
        let mut dna = dna.into_iter();
        let frames = std::array::from_fn(|_| {
            let iter = dna.clone().trimmed_to_codon();
            let _ = dna.next();
            iter
        });

        let ncbi1 = TranslationTable::Ncbi1.to_fn();
        let aas = frames.clone().map(|frame| {
            let translated: Vec<_> = frame.codons().map(ncbi1).collect();
            String::from_utf8(translated).unwrap()
        });
        let aas_rc = frames.map(|frame| {
            let translated: Vec<_> = frame.reverse_complement().codons().map(ncbi1).collect();
            String::from_utf8(translated).unwrap()
        });

        Self { aas, aas_rc }
    }

    /// Returns an iterator over the windows and their original position in the full FASTA
    ///
    /// The indexes returned by `enumerate_windows` are not strictly
    /// monotonically increasing because the current implementation (e.g.)
    /// yields all forward frame 0 windows before all RC frame 0 windows.
    pub fn enumerate_windows(
        &self,
        aa_window_len: usize,
    ) -> impl ExactSizeIterator<Item = (usize, &str)> {
        let forward_aas = |offset: usize| {
            ascii_str_windows(&self.aas[offset], aa_window_len)
                .enumerate()
                .map(move |(i, w)| (offset + 3 * i, w))
        };
        let reverse_aas = |offset: usize| {
            ascii_str_windows(&self.aas_rc[offset], aa_window_len)
                .rev()
                .enumerate()
                .map(move |(i, w)| (offset + 3 * i, w))
        };

        WithExactSize(
            // Deliberately using chain instead of flat_map for accurate size_hints
            forward_aas(0)
                .chain(forward_aas(1))
                .chain(forward_aas(2))
                .chain(reverse_aas(0))
                .chain(reverse_aas(1))
                .chain(reverse_aas(2)),
        )
    }
}

// Unfortunate, but needed as long as we're using strings, AFAICT
fn ascii_str_windows(
    data: &str,
    window_len: usize,
) -> impl DoubleEndedIterator<Item = &str> + ExactSizeIterator {
    let num_windows = data.as_bytes().windows(window_len).len();
    (0..num_windows).map(move |i| &data[i..i + window_len])
}

// Workaround for chains not having exact sizes even if they wouldn't overflow.
struct WithExactSize<I>(I);

impl<I: Iterator> Iterator for WithExactSize<I> {
    type Item = I::Item;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<I: Iterator> ExactSizeIterator for WithExactSize<I> {}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use quickcheck::quickcheck;
    use quickdna::{BaseSequence, DnaSequence, Nucleotide};

    use super::*;

    macro_rules! assert_contains {
        ($set:expr, $value:expr) => {
            let set = &$set;
            let value = &$value;
            assert!(set.contains(value), "{value:?} not in {set:?}");
        };
    }

    fn to_dna(repr: &str) -> Vec<Nucleotide> {
        let dna: DnaSequence<Nucleotide> = repr.parse().unwrap();
        dna.as_slice().to_owned()
    }

    fn expected_dna_windows(dna: &[Nucleotide], window_len: usize) -> HashSet<(usize, String)> {
        let mut expected = HashSet::new();

        for (frame_idx, window) in dna.windows(window_len).enumerate() {
            let window_dna = DnaSequence::<Nucleotide>::new(window.to_vec());
            expected.insert((frame_idx, window_dna.to_string()));
            expected.insert((frame_idx, window_dna.reverse_complement().to_string()));
        }

        expected
    }

    fn expected_protein_windows(dna: &[Nucleotide], window_len: usize) -> HashSet<(usize, String)> {
        let mut expected = HashSet::new();
        let ncbi1 = TranslationTable::Ncbi1;

        for (frame_idx, window) in dna.windows(3 * window_len).enumerate() {
            let window_dna = DnaSequence::<Nucleotide>::new(window.to_vec());
            expected.insert((frame_idx, window_dna.translate(ncbi1).to_string()));
            expected.insert((
                frame_idx,
                window_dna.reverse_complement().translate(ncbi1).to_string(),
            ));
        }

        expected
    }

    #[test]
    fn smoke_test() {
        let dna = to_dna(
            "ATGAGTCTTTTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTC\
             AAAGCCGAGATCGCGCAGAGACTTGAAGATGTCTTTGCAGGGAAGAACACCGATCTTGAG",
        );
        let dna_window_len = 42;
        let protein_window_len = 20;

        let windows = DnaAndProteinWindows::from_dna(dna);

        let vals: HashSet<_> = windows
            .enumerate_windows(dna_window_len, protein_window_len)
            .collect();

        // Windows starting at the first 4 nucleotides

        assert_contains!(vals, (0, "ATGAGTCTTTTAACCGAGGTCGAAACGTACGTTCTCTCTATC")); // DNA
        assert_contains!(vals, (0, "GATAGAGAGAACGTACGTTTCGACCTCGGTTAAAAGACTCAT")); // RC DNA
        assert_contains!(vals, (0, "MSLLTEVETYVLSIVPSGPL")); // AA
        assert_contains!(vals, (0, "EGA*RDDRENVRFDLG*KTH")); // RC AA

        assert_contains!(vals, (1, "TGAGTCTTTTAACCGAGGTCGAAACGTACGTTCTCTCTATCG")); // DNA
        assert_contains!(vals, (1, "CGATAGAGAGAACGTACGTTTCGACCTCGGTTAAAAGACTCA")); // RC DNA
        assert_contains!(vals, (1, "*VF*PRSKRTFSLSSRQAPS")); // AA
        assert_contains!(vals, (1, "*GGLTGR*RERTFRPRLKDS")); // RC AA

        assert_contains!(vals, (2, "GAGTCTTTTAACCGAGGTCGAAACGTACGTTCTCTCTATCGT")); // DNA
        assert_contains!(vals, (2, "ACGATAGAGAGAACGTACGTTTCGACCTCGGTTAAAAGACTC")); // RC DNA
        assert_contains!(vals, (2, "ESFNRGRNVRSLYRPVRPPQ")); // AA
        assert_contains!(vals, (2, "LRGPDGTIERTYVSTSVKRL")); // RC AA

        assert_contains!(vals, (3, "AGTCTTTTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTC")); // DNA
        assert_contains!(vals, (3, "GACGATAGAGAGAACGTACGTTTCGACCTCGGTTAAAAGACT")); // RC DNA
        assert_contains!(vals, (3, "SLLTEVETYVLSIVPSGPLK")); // AA
        assert_contains!(vals, (3, "FEGA*RDDRENVRFDLG*KT")); // RC AA

        // The last 4 AA windows

        assert_contains!(vals, (57, "LKAEIAQRLEDVFAGKNTDL")); // AA
        assert_contains!(vals, (57, "KIGVLPCKDIFKSLRDLGFE")); // RC AA

        assert_contains!(vals, (58, "SKPRSRRDLKMSLQGRTPIL")); // AA
        assert_contains!(vals, (58, "QDRCSSLQRHLQVSARSRL*")); // RC AA

        assert_contains!(vals, (59, "QSRDRAET*RCLCREEHRS*")); // AA
        assert_contains!(vals, (59, "SRSVFFPAKTSSSLCAISAL")); // RC AA

        assert_contains!(vals, (60, "KAEIAQRLEDVFAGKNTDLE")); // AA
        assert_contains!(vals, (60, "LKIGVLPCKDIFKSLRDLGF")); // RC AA

        // The last 4 DNA windows

        assert_contains!(vals, (75, "CAGAGACTTGAAGATGTCTTTGCAGGGAAGAACACCGATCTT")); // DNA
        assert_contains!(vals, (75, "AAGATCGGTGTTCTTCCCTGCAAAGACATCTTCAAGTCTCTG")); // RC DNA

        assert_contains!(vals, (76, "AGAGACTTGAAGATGTCTTTGCAGGGAAGAACACCGATCTTG")); // DNA
        assert_contains!(vals, (76, "CAAGATCGGTGTTCTTCCCTGCAAAGACATCTTCAAGTCTCT")); // RC DNA

        assert_contains!(vals, (77, "GAGACTTGAAGATGTCTTTGCAGGGAAGAACACCGATCTTGA")); // DNA
        assert_contains!(vals, (77, "TCAAGATCGGTGTTCTTCCCTGCAAAGACATCTTCAAGTCTC")); // RC DNA

        assert_contains!(vals, (78, "AGACTTGAAGATGTCTTTGCAGGGAAGAACACCGATCTTGAG")); // DNA
        assert_contains!(vals, (78, "CTCAAGATCGGTGTTCTTCCCTGCAAAGACATCTTCAAGTCT")); // RC DNA

        assert_eq!(
            vals.len(),
            // The full DNA above should have a length of 120.
            // For DNA windows: after the first 42 bp window there should be 78 more bps, each of
            // which allows another window to exist, for a total of 79 windows.
            79
            // For protein windows: A protein window is 20 aa or 60 bp, after which there are
            // 60 more bp, for a total 61 windows.
            + 61
            // RC doesn't change DNA length, so the number of RC windows should be the same.
            + 79 + 61
        );
    }

    #[test]
    fn expected_windows() {
        let dna = to_dna(
            "ATGAGTCTTTTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTC\
             AAAGCCGAGATCGCGCAGAGACTTGAAGATGTCTTTGCAGGGAAGAACACCGATCTTGAG",
        );
        let dna_window_len = 42;
        let protein_window_len = 20;

        let windows = DnaAndProteinWindows::from_dna(&dna);
        let enumerated = windows.enumerate_windows(dna_window_len, protein_window_len);

        assert_eq!(
            enumerated.len(),
            2 * (dna.len() - dna_window_len + 1) + 2 * (dna.len() - 3 * protein_window_len + 1)
        );

        let actual = enumerated
            .map(|(idx, value)| (idx, value.to_owned()))
            .collect::<HashSet<(usize, String)>>();

        let expected = &expected_dna_windows(&dna, dna_window_len)
            | &expected_protein_windows(&dna, protein_window_len);

        // Note: assert_eq! produces overwhelming failure messages; this shows helpful diffs
        let missing = &expected - &actual;
        assert!(missing.is_empty(), "missing (index, window)s: {missing:?}");
        let unexpected = &actual - &expected;
        assert!(
            unexpected.is_empty(),
            "unexpected (index, window)s: {unexpected:?}"
        );
    }

    #[test]
    fn windows_convenience_functions() {
        let dna_window_len = 42;
        let protein_window_len = 20;

        {
            let dna = to_dna(
                "ATGAGTCTTTTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTC\
                 AAAGCCGAGATCGCGCAGAGACTTGAAGATGTCTTTGCAGGGAAGAACACCGATCTTGAG",
            );

            let windows = DnaAndProteinWindows::from_dna(&dna);
            let enumerated = windows.enumerate_windows(dna_window_len, protein_window_len);
            let protein_windows = ProteinWindows::from_dna(&dna);
            let enumerated_protein_windows = protein_windows.enumerate_windows(protein_window_len);

            assert_eq!(enumerated_protein_windows.len(), 122);
            assert_eq!(enumerated.len(), 280);
        }

        {
            let dna = to_dna("ATG");

            let windows = DnaAndProteinWindows::from_dna(&dna);
            let enumerated = windows.enumerate_windows(dna_window_len, protein_window_len);
            let protein_windows = ProteinWindows::from_dna(&dna);
            let enumerated_protein_windows = protein_windows.enumerate_windows(protein_window_len);

            assert_eq!(enumerated_protein_windows.len(), 0);
            assert_eq!(enumerated.len(), 0);
        }

        {
            let dna = to_dna("ATGAGTCTTTTAACCGAGGTCGAAACGTACGTTCTCTCTATC");

            let windows = DnaAndProteinWindows::from_dna(dna);
            let enumerated = windows.enumerate_windows(dna_window_len, protein_window_len);

            assert_eq!(enumerated.len(), 2);
        }
    }

    #[test]
    #[should_panic]
    fn test_dna_and_protein_windows_cannot_have_zero_dna_len() {
        let dna = to_dna("CATTAGATCG");
        let dna_window_len = 0;
        let protein_window_len = 3;
        let windows = DnaAndProteinWindows::from_dna(dna);
        windows
            .enumerate_windows(dna_window_len, protein_window_len)
            .count();
    }

    #[test]
    #[should_panic]
    fn test_dna_and_protein_windows_cannot_have_zero_protein_len() {
        let dna = to_dna("CATTAGATCG");
        let dna_window_len = 6;
        let protein_window_len = 0;
        let windows = DnaAndProteinWindows::from_dna(dna);
        windows
            .enumerate_windows(dna_window_len, protein_window_len)
            .count();
    }

    #[test]
    #[should_panic]
    fn test_dna_windows_cannot_have_zero_len() {
        let dna = to_dna("CATTAG");
        let window_len = 0;
        let windows = DnaWindows::from_dna(dna.as_slice());
        windows.enumerate_windows(window_len).count();
    }

    #[test]
    #[should_panic]
    fn test_protein_windows_cannot_have_zero_len() {
        let dna = to_dna("CATTAGATCG");
        let window_len = 0;
        let windows = ProteinWindows::from_dna(dna);
        windows.enumerate_windows(window_len).count();
    }

    #[test]
    fn test_dna_and_protein_windows_len_can_be_calculated() {
        let dna = to_dna("CATTAGATCG");
        let dna_window_len = 9;
        let protein_window_len = 2;
        let windows = DnaAndProteinWindows::from_dna(dna);

        // Doing this to make sure the WithExactSize wrapper works.
        let len = windows
            .enumerate_windows(dna_window_len, protein_window_len)
            .len();
        assert_eq!(len, 14); // 2 forward dna, 2 rc dna, 5 forward protein, 5 rc protein
    }

    #[test]
    fn test_dna_windows_len_can_be_calculated() {
        let dna = to_dna("CATTAG");
        let window_len = 3;
        let windows = DnaWindows::from_dna(dna);

        // Doing this to make sure the WithExactSize wrapper works.
        let len = windows.enumerate_windows(window_len).len();
        assert_eq!(len, 8); // 4 forward, 4 rc
    }

    #[test]
    fn test_protein_windows_len_can_be_calculated() {
        let dna = to_dna("CATTAGATCG");
        let window_len = 2;
        let windows = ProteinWindows::from_dna(dna);

        // Doing this to make sure the WithExactSize wrapper works.
        let len = windows.enumerate_windows(window_len).len();
        assert_eq!(len, 10); // 5 forward, 5 rc
    }

    quickcheck! {

        fn dna_and_protein_windows_matches_reference_implementation(dna: Vec<Nucleotide>) -> bool {
            let dna_window_len = 42;
            let protein_window_len = 20;
            let windows = DnaAndProteinWindows::from_dna(&dna);

            let actual: HashSet<_> = windows
                .enumerate_windows(dna_window_len, protein_window_len)
                .map(|(idx, value)| (idx, value.to_owned()))
                .collect();

            let expected = &expected_dna_windows(&dna, dna_window_len)
                | &expected_protein_windows(&dna, protein_window_len);

            actual == expected
        }

        fn dna_windows_matches_reference_implementation(dna: Vec<Nucleotide>) -> bool {
            let window_len = 42;
            let windows = DnaWindows::from_dna(&dna);

            let actual: HashSet<_> = windows
                .enumerate_windows(window_len)
                .map(|(idx, value)| (idx, value.to_owned()))
                .collect();

            let expected = expected_dna_windows(&dna, window_len);

            actual == expected
        }

        fn protein_windows_matches_reference_implementation(dna: Vec<Nucleotide>) -> bool {
            let window_len = 20;
            let windows = ProteinWindows::from_dna(&dna);

            let actual: HashSet<_> = windows
                .enumerate_windows(window_len)
                .map(|(idx, value)| (idx, value.to_owned()))
                .collect();

            let expected = expected_protein_windows(&dna, window_len);

            actual == expected
        }
    }
}
