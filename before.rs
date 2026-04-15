use quickdna::{BaseSequence, DnaSequence, FastaRecord, NucleotideLike, TranslationTable};

use crate::parsefasta::{DNA_WINDOW_LEN, PROTEIN_WINDOW_LEN};

#[derive(Debug, Clone)]
pub struct SequenceWindows {
    pub dna_sequence_length: usize,
    pub dna: Vec<String>,
    pub dna_reverse: Vec<String>,
    pub aa: Vec<Vec<String>>,
}

impl SequenceWindows {
    /// returns the combined length of all types of windows
    pub fn len(&self) -> usize {
        self.dna.len()
            + self.dna_reverse.len()
            + self
                .aa
                .iter()
                .map(|ele| ele.len())
                .reduce(|a, b| a + b)
                .unwrap_or(0)
    }

    /// returns 'true' if all window vectors are empty
    pub fn is_empty(&self) -> bool {
        self.dna.is_empty() && self.dna_reverse.is_empty() && self.aa.iter().all(|a| a.is_empty())
    }

    /// Returns an iterator over the windows and their original position in the full FASTA
    /// The indexes returned by enumerate are not strictly monotonically increasing
    pub fn enumerate(&self) -> impl Iterator<Item = (usize, &str)> {
        self.dna
            .iter()
            .enumerate()
            .map(|(i, s)| (i, s.as_str()))
            .chain(
                self.dna_reverse
                    .iter()
                    .enumerate()
                    .map(|(index, value)| (self.dna_reverse.len() - index - 1, value.as_str())),
            )
            .chain(
                self.aa
                    .iter()
                    .enumerate()
                    .flat_map(move |(operation, frame_windows)| {
                        // this code is the reverse of QuickDNA translate_all_frames
                        // https://github.com/SecureDNA/quickdna/blob/main/src/rust_api.rs#L196
                        let dna_per_aa = 3;
                        frame_windows
                            .iter()
                            .enumerate()
                            .map(move |(index_within_frame, window)| {
                                if operation == 0 || operation == 1 || operation == 2 {
                                    // forward
                                    (dna_per_aa * index_within_frame + operation, window.as_str())
                                } else {
                                    // backwards
                                    let index = self.dna_sequence_length
                                        - ((index_within_frame + PROTEIN_WINDOW_LEN) * dna_per_aa
                                            + (operation % 3));
                                    (index, window.as_str())
                                }
                            })
                    }),
            )
    }

    /// Construct all windows from a FASTA record
    pub fn from_record<T: NucleotideLike>(record: &FastaRecord<DnaSequence<T>>) -> SequenceWindows {
        let translations = record
            .contents
            .translate_all_frames(TranslationTable::Ncbi1);

        let windows = SequenceWindows {
            dna_sequence_length: record.contents.len(),
            dna: record
                .contents
                .windows(DNA_WINDOW_LEN)
                .map(|s| s.to_string())
                .collect(),
            dna_reverse: record
                .contents
                .reverse_complement()
                .windows(DNA_WINDOW_LEN)
                .map(|s| s.to_string())
                .collect(),
            aa: translations
                .iter()
                .map(|t| {
                    t.windows(PROTEIN_WINDOW_LEN)
                        .map(|s| s.to_string())
                        .collect()
                })
                .collect(),
        };

        // amino acids
        windows
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::str::FromStr;

    use quickdna::{FastaParseSettings, FastaParser, Nucleotide};

    use super::*;

    #[test]
    fn expected_windows() {
        let order_fasta = "ATGAGTCTTTTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCGCAGAGACTTGAAGATGTCTTTGCAGGGAAGAACACCGATCTTGAG";

        let parser = FastaParser::<DnaSequence<Nucleotide>>::new(
            FastaParseSettings::new()
                .concatenate_headers(true)
                .allow_preceding_comment(false),
        );

        let fastas = parser.parse_str(order_fasta).unwrap();

        let windows = SequenceWindows::from_record(&fastas.records[0]);

        assert_eq!(windows.dna.len(), order_fasta.len() - DNA_WINDOW_LEN + 1);
        assert_eq!(
            windows.dna_reverse.len(),
            order_fasta.len() - DNA_WINDOW_LEN + 1
        );

        let m = windows
            .enumerate()
            .map(|(idx, value)| (value.to_owned(), idx))
            .collect::<HashMap<String, usize>>();

        for frame_idx in 0..order_fasta.len() - PROTEIN_WINDOW_LEN * 3 {
            let translation = DnaSequence::<Nucleotide>::from_str(
                &order_fasta[frame_idx..frame_idx + PROTEIN_WINDOW_LEN * 3],
            )
            .unwrap()
            .translate(TranslationTable::Ncbi1);
            assert_eq!(m[&translation.to_string()], frame_idx)
        }

        for frame_idx in 0..order_fasta.len() - PROTEIN_WINDOW_LEN * 3 {
            let translation = DnaSequence::<Nucleotide>::from_str(
                &order_fasta[frame_idx..frame_idx + PROTEIN_WINDOW_LEN * 3],
            )
            .unwrap()
            .reverse_complement()
            .translate(TranslationTable::Ncbi1);
            assert_eq!(
                m[&translation.to_string()],
                frame_idx,
                "{} was not where expected",
                translation
            )
        }

        for frame_idx in 0..order_fasta.len() - DNA_WINDOW_LEN {
            let translation = DnaSequence::<Nucleotide>::from_str(
                &order_fasta[frame_idx..frame_idx + DNA_WINDOW_LEN],
            )
            .unwrap();
            assert_eq!(m[&translation.to_string()], frame_idx)
        }

        for frame_idx in 0..order_fasta.len() - DNA_WINDOW_LEN {
            let translation = DnaSequence::<Nucleotide>::from_str(
                &order_fasta[frame_idx..frame_idx + DNA_WINDOW_LEN],
            )
            .unwrap()
            .reverse_complement();
            assert_eq!(m[&translation.to_string()], frame_idx)
        }
    }

    #[test]
    fn windows_convenience_functions() {
        {
            let order_fasta = "ATGAGTCTTTTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCGCAGAGACTTGAAGATGTCTTTGCAGGGAAGAACACCGATCTTGAG";

            let parser = FastaParser::<DnaSequence<Nucleotide>>::new(
                FastaParseSettings::new()
                    .concatenate_headers(true)
                    .allow_preceding_comment(false),
            );

            let fastas = parser.parse_str(order_fasta).unwrap();

            let windows = SequenceWindows::from_record(&fastas.records[0]);
            assert_eq!(windows.aa.len(), 6);
            assert!(!windows.is_empty());
            assert_eq!(windows.len(), 280);
        }

        {
            let order_fasta = "ATG";

            let parser = FastaParser::<DnaSequence<Nucleotide>>::new(
                FastaParseSettings::new()
                    .concatenate_headers(true)
                    .allow_preceding_comment(false),
            );

            let fastas = parser.parse_str(order_fasta).unwrap();

            let windows = SequenceWindows::from_record(&fastas.records[0]);

            assert_eq!(windows.aa.len(), 2);
            assert!(windows.is_empty());
            assert_eq!(windows.len(), 0);
        }

        {
            let order_fasta = "ATGAGTCTTTTAACCGAGGTCGAAACGTACGTTCTCTCTATC";

            let parser = FastaParser::<DnaSequence<Nucleotide>>::new(
                FastaParseSettings::new()
                    .concatenate_headers(true)
                    .allow_preceding_comment(false),
            );

            let fastas = parser.parse_str(order_fasta).unwrap();

            let windows = SequenceWindows::from_record(&fastas.records[0]);
            assert!(!windows.is_empty());
            assert_eq!(windows.len(), 2);
        }
    }
}
