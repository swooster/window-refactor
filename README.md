# Example refactor of SecureDNA windowing code

This repo shows a refactor I'm particularly proud of because it simultaneously
made the code easier to reason about while reducing memory use.

[(before)](./before.rs)

[(after)](./after.rs)


## Purpose of code

At a high level, SecureDNA's software scans DNA for hazards by splitting it
into overlapping windows, hashing and encrypting the windows and comparing
the hashes against pre-encrypted hazard hashes. An added wrinkle is that
there are various transformations to the DNA that must be checked, such as
reverse-complementation and translation of DNA into peptides.

The snippets of code in this repo were used to build the various
(transformed) windows before they're hashed.


## Tiny biology primer

(This skips describing RNA because it's not necessary to understand
the code.)

DNA is a sequence of nucleotides. A nucleotide can be any one of
`A` (adenine), `C` (cytosine), `G` (guanine) or `T` (thymine).

You can "complement" DNA by swapping `A`s with `T`s (and vice versa)
as well as `C`s with `G`s (and vice versa). For example, the complement
of `GAGACATTAG` is `CTCTGTAATC`. A very common operation for DNA is
the "reverse complement", which consists of... reversing and complementing
it. The reverse complement of `GAGACATTAG` is `CTAATGTCTC`.

DNA can be translated into sequences of amino acids (called "peptides" or
"proteins"). To do so, you chunk it into nucleotide triplets (called "codons")
and apply a mapping from codons to amino acids. The most common mapping is
[NCBI1](https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables#Standard_RNA_codon_table).
Whereas there are only 4 kinds of nucleotides, there are about 20-22 amino
acids (each with their own letter abbreviation), plus "stop" (written as `*`).
Crucially, because you're chunking into triplets, there are three possible
starting offsets (called "reading frames") that produce very different
peptides; a starting offset of 3 would produce a peptide that's a suffix of
the peptide produced by a starting offset of 0, so they're considered
to be the same reading frame. For example the three reading frames of
`GAGACATTAG` produce `ETL`, `RH*` and `DI` when mapped through NCBI1 because:
* `GAG ACA TTA G` &rarr; `GAG ACA TTA` &rarr; `ETL`
* `G AGA CAT TAG` &rarr; ` AGA CAT TAG` &rarr; `RH*`
* `GA GAC ATT AG` &rarr; `  GAC ATT` &rarr; `DI`


## What the code actually does

The goal of the code is pretty easy: Given a DNA sequence `dna`, and lengths `D`
and `P`, the code generates the following data: (though not in this order)
* For `start` in `0..dna.len()`:
  * If `&dna[start..start+D]` is a valid slice:
    * Yield `(start, dna[start..start+D])`
    * Yield `(start, reverse_complement(dna[start..start+D]))`
  * If `&dna[start..start+P]` is a valid slice:
    * Yield `(start, ncbi1(dna[start..start+P]))`
    * Yield `(start, ncbi1(reverse_complement(dna[start..start+P])))`

(all windows are encoded via their ASCII representation)

The complexity comes from trying to deduplicate the work of applying
transformations like reverse-complement or mapping though NCBI1.

## Issues with the previous version

### Performance

The previous code
[generates all windows preemptively](https://github.com/swooster/window-refactor/blob/main/before.rs#L73-L98)
without using references. This means that memory usage is
`O(length of DNA × (length of DNA windows + length of peptide windows))`.
I'll spare you the gory details, but with a DNA window length of 42
nucleotides and peptide window length of 60 nucleotides (or 20 amino
acids), a 1MB DNA sequence can easily balloon to over 160MB of window
data, with an additional 40MB wasted on unused capacity.

Another thing worth pointing out is that giving each window its own
allocation is bad for memory locality. Even if the memory allocator
happened to place windows contiguously, the duplication of overlapping
data still puts unnecessary pressure on the cache.


### Reasoning


The code needs to yield many `(start, window)` items and
the calculations for `start` are hard to verify. They're
[initially trivial](https://github.com/swooster/window-refactor/blob/main/before.rs#L36-L37)
for forward DNA windows, then
[require attention](https://github.com/swooster/window-refactor/blob/main/before.rs#L41-L42)
for reverse complement DNA windows (Did you catch that
`self.dna_reverse.len()` isn't the same as `self.dna_sequence_length`?
Would you have noticed if the `- 1` were  a `+ 1`?), then
[become exceedingly tricky to get right](https://github.com/swooster/window-refactor/blob/main/before.rs#L47-L67)
for the peptide windows. In fact, whether the code is correct
depends on undocumented details of how
[`quickdna::DnaSequence::translate_all_frames`](https://github.com/SecureDNA/quickdna/blob/ad2129e9e0e4bbbc08caf520cfbbc9bd11682d43/src/rust_api.rs#L262-L275)
defines reverse complement reading frames (do you measure the
offsets from the beginning of the reverse complemented sequence
or the original sequence?). The fiddly details of the code made the
original author quite hesitant about the prospect of letting anyone
refactor it.

## The refactor

There are a couple of things I'd handle differently if I had a do-over
(`WithExactSize` and the overreliance on iterators for input).
I'll come back to those later.

### Performance

The refactored code avoids precomputing windows. Rather, it stores
the original DNA, the reverse complement of the DNA, all 3 forward
translations of the DNA and all 3 reverse complement translations
of the DNA, so memory usage is just `O(length of DNA)`. A 1MB DNA
sequence inflates to 4MB when all transformations are taken into
account. Yielded windows are just (slightly shifting) slices of
the backing buffers so they don't cost much extra memory and
exhibit good locality.

### Reasoning

This is the biggest difference. I feel that the refactored code
is a lot easier to verify.

Starting at a low level there's
[`ascii_str_windows`](https://github.com/swooster/window-refactor/blob/main/after.rs#L153-L160),
which is like
[`[T]::windows`](https://doc.rust-lang.org/std/primitive.slice.html#method.windows)
for `str`s with ASCII-only content. One potentially odd choice
is that I eschewed `data.len() - window_len + 1` in favor of
`data.as_bytes().windows(window_len).len()`. The reason is that,
whether used for calculating indices or the number of windows,
`L - W + 1` is surprisingly error-prone for such a simple
expression; both the orginal code and the in-development refactored
code encountered off-by-one bugs stemming from it. In contrast,
`num_windows` is `data.as_bytes().windows(window_len).len()`
by definition, and
[that compiles down to similar machine code](https://godbolt.org/z/8hPGG748v)
(except that it asserts `window_len != 0`).

Moving up a level to
[DnaWindows::enumerate_windows](https://github.com/swooster/window-refactor/blob/main/after.rs#L69-L83),
the window `start` calculation is
[trivial for forward windows](https://github.com/swooster/window-refactor/blob/main/after.rs#L78) and
[easy for reverse-complement windows](https://github.com/swooster/window-refactor/blob/main/after.rs#L79-L81):
because reverse complement reverses the DNA,
`ascii_str_windows(&self.dna_rc, window_len)` yields windows from back
to front (relative to the original DNA sequence) and reversing the
window order is all it takes for `enumerate()` to do the right thing.

Protein window calculations are a little more complex, but much simpler
than before. The
[constructor for `ProteinWindows`](https://github.com/swooster/window-refactor/blob/main/before.rs#L71-L102)
avoids relying on 
[`quickdna::DnaSequence::translate_all_frames`](https://github.com/SecureDNA/quickdna/blob/ad2129e9e0e4bbbc08caf520cfbbc9bd11682d43/src/rust_api.rs#L262-L275)
due to that methods's tricky behavior around starting offsets for
reverse-complement reading frames (they depend on the length of the DNA).
It instead
[populates `aas` and `aas_rc`](https://github.com/swooster/window-refactor/blob/main/after.rs#L108-L115)
(short for "amino acids" and "amino acids, reverse complemented")
such that source DNA for `aas_rc[reading_frame]` is the reverse
complement of the source DNA for `aas[reading_frame]`, and they
both have a starting offset of `reading_frame` in the original
DNA sequence.

The [logic for the (forward) peptides](https://github.com/swooster/window-refactor/blob/main/after.rs#L129-L133)
isn't too bad; the big tweak is the window `start` for `aas[reading_frame]`
begins at `reading_frame` and increases by 3 nucleotides every time the peptide
window shifts by an amino acid, hence
[this index adjustment](https://github.com/swooster/window-refactor/blob/main/after.rs#L132).
The [logic for reverse complement peptides](https://github.com/swooster/window-refactor/blob/main/after.rs#L134-L139)
is the same idea, except that
[the window order needs to be reversed](https://github.com/swooster/window-refactor/blob/main/after.rs#L136)
for the exact same reason that
[the reverse complement DNA window order needs to be reversed](https://github.com/swooster/window-refactor/blob/main/after.rs#L80).

Another improvement worth noting is how the total number of windows
are calculated. The original file has separate code for
[enumerating windows](https://github.com/swooster/window-refactor/blob/main/before.rs#L31-L69)
and [counting windows](https://github.com/swooster/window-refactor/blob/main/before.rs#L14-L24).
The refactored file doesn't have separate code that can fall out of
sync; it just relies on the standard library's iterator adaptors
knowing how to properly measure their own length, so the calculated
total number of windows inherently matches the definition of the iterator.

### Downsides of refactor

There were a couple of things I think I got wrong with this refactor.

First and foremost,
[`WithExactSize`](https://github.com/swooster/window-refactor/blob/main/after.rs#L162-L178)
can in theory panic; I had figured memory limitations elsewhere would
prevent that, but reducing the need for non-local reasoning is a big
part of what makes Rust so delightful, so I think this was a bad choice.
The silver lining is that it won't silently overflow.

The other place where I think I mistepped was overly relying on
iterators for the sequence input. In hindsight, the system didn't
need that much generality, so it was an unnecessary burden on the
readablility of the code.

## Summary

With the refactor, it's easier to verify that the window start indices
are correct, initial setup takes less processing, the window
generation code holds 50x less memory, and iterating over windows is
nicer on the cache.

## Copyright

The rust code in this repository is copyrighted by
[SecureDNA Stiftung (SecureDNA Foundation)](mailto:licensing@securedna.org)
and shared from the private repo with permission:
* [`before.rs`](./before.rs) is from `47160640368ec90c227eb7161604c99b4f967d8e:synth_client/src/windows.rs` (2023-01-06).
* [`after.rs`](./after.rs) is from `07437d74f9a4ea9212406ebdff89d2ddd2f581cc:synth_client/src/windows.rs` (2023-01-30).

The public releases repo can be found at: <https://github.com/SecureDNA/SecureDNA>
