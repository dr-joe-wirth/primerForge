<img src="assets/primerforge.png" alt="primerForge" width="250" height="250">

[![DOI](https://joss.theoj.org/papers/10.21105/joss.06850/status.svg)](https://doi.org/10.21105/joss.06850)

# primerForge

software to identify primers that can be used to distinguish genomes

If you use this software, please cite [our article](https://doi.org/10.21105/joss.06850) in the Journal of Open Source Software.

## Installation
> [!NOTE]
> `primerForge` is incompatible with versions python3.8 and below and python3.12 and above.

### `pip` installation
`primerForge` can be installed with the following commands:

```shell
pip install primerforge
conda install ispcr
```

### `pixi` installation
> [!NOTE]
> Pixi is experimental

To install with Pixi, download [`pixi.toml`](pixi.toml) and then run the following commands:

```bash
pixi shell
```


### `conda` installation
A `conda` installation is currently unavailable due to `khmer` being unsupported. We are actively working to resolve this.

### Manual installation

> [!NOTE]
> This might take up to ten minutes.

```shell
git clone https://github.com/dr-joe-wirth/primerForge.git
conda env create -f primerForge/environment.yml
conda activate primerforge
```


### Docker installation

A Docker image for the latest release is available at [DockerHub](https://hub.docker.com/r/jwirth/primerforge)

### Checking installation

If `primerForge` is installed correctly, then the following command should execute without errors:

```shell
primerForge --check_install
```

If you installed manually, you may need to use the following command instead:

```shell
python primerForge.py --check_install
```

### Running unit tests
In order to run unit tests, install `primerForge` using [the instructions above](#installation). You will also need to clone the repository if you haven't already:

```bash
git clone https://github.com/dr-joe-wirth/primerForge.git
```

Once installed and cloned, run the following commands to run the unit tests:

```bash
python3 -m unittest discover -s ./primerForge/bin/unit_tests/ -p "*_test.py"
```


## Usage
### Basic usage

```text
usage:
    primerForge [-ioBubfpgtrdnkvh]

required arguments:
    -i, --ingroup          [file] ingroup filename or a file pattern inside double-quotes (eg."*.gbff")

optional arguments:
    -o, --out              [file] output filename for primer pair data (default: results.tsv)
    -B, --bed_file         [file] output filename for primer data in BED file format (default: primers.bed)
    -u, --outgroup         [file] outgroup filename or a file pattern inside double-quotes (eg."*.gbff")
    -b, --bad_sizes        [int,int] a range of PCR product lengths that the outgroup cannot produce (default: same as '--pcr_prod')
    -f, --format           [str] file format of the ingroup and outgroup [genbank|fasta] (default: genbank)
    -p, --primer_len       [int(s)] a single primer length or a range specified as 'min,max'; (minimum 10) (default: 16,20)
    -g, --gc_range         [float,float] a min and max percent GC specified as a comma separated list (default: 40.0,60.0)
    -t, --tm_range         [float,float] a min and max melting temp (°C) specified as a comma separated list (default: 55.0,68.0)
    -r, --pcr_prod         [int(s)] a single PCR product length or a range specified as 'min,max' (default: 120,2400)
    -d, --tm_diff          [float] the maximum allowable Tm difference °C between a pair of primers (default: 5.0)
    -n, --num_threads      [int] the number of threads for parallel processing (default: 1)
    -k, --keep             keep intermediate files (default: False)
    -v, --version          print the version
    -h, --help             print this message

    --check_install        check installation
    --debug                run in debug mode
    --advanced             print advanced options
```

### Advanced options
#### `primer3` parameters
```text
--primer3_mv_conc      [float] monovalent cation concentration (mM) (default: 50.0)
--primer3_dv_conc      [float] divalent cation concentration (mM) (default: 1.5)
--primer3_dntp_conc    [float] dNTP concentration (mM) (default: 0.6)
--primer3_dna_conc     [float] template DNA concentration (nM) (default: 50.0)
--primer3_temp_c       [float] simulation temp (°C) for ΔG (default: 37.0)
--primer3_max_loop     [int] maximum size (bp) of loops in primer secondary structures (default: 30)
```

#### `isPcr` parameters
```text
--isPcr_minGood        [int] minimum size (bp) where there must be 2 matches for each mismatch (default: 6)
--isPcr_minPerfect     [int] minimum size (bp) of perfect match at 3' end of primer (default: 8)
--isPcr_tileSize       [int] the size of match that triggers an alignment (default: 10)
```

#### Additional `primerForge` parameters
```text
--temp_tolerance       [float] minimum number of degrees (°C) below primer Tm allowed for secondary structure Tm (default: 5.0)
--max_repeats          [int] maximum allowed length (bp) of homopolymers (repeats) in primer sequences (default: 4)
--bin_size             [int] maximum allowed length (bp) of contiguous regions of overlapping primers (bins) (default: 64)
```

## Results
### `results.tsv`
The results produced by `primerForge` list one primer pair per line, and these pairs are sorted based on the following criteria (in order of importance):

  * largest difference from the ingroup PCR product range for outgroup PCR products
  * lowest number of outgroup PCR products
  * lowest G+C difference between the pair
  * lowest Tm difference between the pair
  * lowest Tm for heterodimers
  * lowest Tm (sum) for homodimers
  * lowest Tm (sum) for hairpins
  * lowest variance in ingroup PCR product sizes
  * largest median ingroup PCR product size

See [the example](#examining-the-results) for details about the file format.

### `primers.bed`
The BED file produced by `primerForge` has the following format:

|column number|title|definition|
|:-----------:|:---:|:----------------|
|1|contig|the name of the contig|
|2|start|the start coordinate of the primer|
|3|end|the end coordinate of the primer (_exclusive_)|
|4|sequence|the sequence of the primer|
|5|pair number|an integer that links forward and reverse primers; overloading the "score" field traditionally found in BED file format|
|6|strand|the DNA strand ('+' or '-')|


## Workflow

```mermaid
flowchart TB
    ingroup[/"ingroup genomes"/]
    ingroup --> A

    %% get unique kmers
    subgraph A["for each genome"]
        uniqKmer["get unique kmers"]
    end

    %% get shared kmers
    sharedKmers(["shared kmers"])
    uniqKmer -- intersection --> sharedKmers

    %% get candidate kmers
    subgraph B["for each genome"]
        subgraph B0["for each kmer start position"]
            subgraph B1["pick one kmer"]
                GC{"GC in
                 range?"}
                Tm{"Tm in
                range?"}
                homo{"repeats
                ≤ 3bp?"}
                hair{"no hairpins?"}
                dime{"no homo-
                dimers?"}
                GC-->Tm-->homo-->hair-->dime
            end
        end
    end

    %% connections up to candidate kmers
    sharedKmers --> B
    dump1[/"checkpoint"/]
    sharedKmers --> dump1
    candidates(["unique, shared kmers; one per start position"])
    dime --> candidates

    %% get primer pairs
    subgraph C["for one genome"]
        bin1["bin overlapping kmers (64bp max)"]
        bin2["remove kmers that are
        substrings of other kmers"]
        bin3["get bin pairs"]

        %% evaluate one kmer pair
        subgraph C0["for each bin pair"]
            size{"is PCR
            size ok?"}
            subgraph C1["for each primer pair"]
                prime{"is 3' end
                G or C?"}
                temp{"is Tm
                difference ok?"}
                gc2{"is GC
                difference ok?"}
                hetero{"no hetero-
                dimers?"}
            end
            size --> C1
        end
    end

    candPair(["candidate primer pairs"])
    allSharePair(["all shared primer pairs"])

    %% get shared primer pairs
    subgraph D["for each candidate primer pair"]
        subgraph D0["for each other genome"]
            pcr{"is PCR
            size ok?"}
        end
    end

    bin1 --> bin2
    bin2 --> bin3
    bin3 --> C0
    prime --> temp --> gc2 --> hetero --> candPair
    candPair --> D
    pcr --> allSharePair

    %% one pair per bin pair
    subgraph E["for each bin pair"]
        keep["keep only one primer pair"]
    end
    
    selectedSharePair(["selected shared
                        primer pairs"])
    dump2[/"checkpoint"/]
    dump3[/"checkpoint"/]

    candidates --> dump2
    candidates --> C
    allSharePair --> E
    keep --> selectedSharePair
    selectedSharePair --> dump3

    %% outgroup removal
    outgroup[/"outgroup genomes"/]

    subgraph F["for each outgroup genome"]
        subgraph F0["for each primer pair"]
            ogsize{"PCR size outside
            disallowed range?"}
        end
    end

    selectedSharePair --> F0
    outgroup --> F
    
    noout(["primer pairs absent from outgroup"])
    ogsize --> noout
    noout --> dump4
    dump4[/"checkpoint"/]

    ispcr["filter pairs using isPcr"]
    final(["final primer pairs"])
    dump5[/"checkpoint"/]
    write[/"sort pairs and
            write to files"/]

    noout --> ispcr
    ispcr --> final
    final --> dump5
    final --> write
```

## Example using _Mycloplasma mycoides_ genomes
This example assumes you have already installed `primerForge` [as described above](#installation).

### Motivation
In this example, we will use `primerForge` to find pairs of primers between 18bp and 24bp that can be used to differentiate three strains of _Mycoplasma mycoides_ subspecies mycoides (the "ingroup") from two strains of _Mycoplasma mycoides_ subspecies capri (the "outgroup"). The primer pairs identified by `primerForge` are predicted to produce a single PCR product between 500bp and 2000bp in the ingroup. These same primer pairs are predicted to produce PCR products <300bp, >3000bp, or no PCR products at all in the outgroup.

### Preparing the workspace
In order to get started, create a directory called `mycoplasma_test` and move into it:

```bash
mkdir ./mycoplasma_test
cd ./mycoplasma_test
```

Next, download the following _Mycoplasma mycoides_ genomes using the following commands:
```bash
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/034/305/GCF_003034305.1_ASM303430v1/GCF_003034305.1_ASM303430v1_genomic.gbff.gz | gzip -d > ./i1.gbff
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/034/275/GCF_003034275.1_ASM303427v1/GCF_003034275.1_ASM303427v1_genomic.gbff.gz | gzip -d > ./i2.gbff
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/034/345/GCF_003034345.1_ASM303434v1/GCF_003034345.1_ASM303434v1_genomic.gbff.gz | gzip -d > ./i3.gbff
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/489/555/GCF_900489555.1_MMC68/GCF_900489555.1_MMC68_genomic.gbff.gz | gzip -d > ./o1.gbff
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/389/745/GCF_018389745.1_ASM1838974v1/GCF_018389745.1_ASM1838974v1_genomic.gbff.gz | gzip -d > ./o2.gbff
```

If you cannot download the genbank files using `curl`, you can download them manually from NCBI by replacing `ftp://` with `http://` and copying and pasting each address into your web browser (eg. `http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/034/305/GCF_003034305.1_ASM303430v1/GCF_003034305.1_ASM303430v1_genomic.gbff.gz`) and then using `gzip -d` on the downloaded file to uncompress it. Finally, be sure to rename each file as shown above (eg. `mv GCF_003034305.1_ASM303430v1_genomic.gbff i1.gbff`).

### Running `primerForge`
We will use the following flags to specify specific parameters for this example:

  * The `--ingroup` and `--outgroup` flags are both file patterns for the ingroup and outgroup genomes, respectively. It is important that this pattern is enclosed in double-quotes as shown below.
  * The `--pcr_prod` flag indicates what sizes we want for ingroup products (500bp to 2,000bp)
  * The `--bad_sizes` flag indicates what sizes we do not want for outgroup products (300bp to 3,000bp).
  * The `--primer_len` flag indicates what length our primers can be (18bp to 24bp)

You can get a list of all available flags using the command `primerForge --help`.

Run `primerForge` using the following command (requires at least 2 Gb of RAM):

```bash
primerForge --ingroup "./i*gbff" --outgroup "./o*gbff" --pcr_prod 500,2000 --bad_sizes 300,3000 --primer_len 18,24
```

After running the command, you should see something like this printed to the screen:

```text
dumping Parameters to 'primerforge_46f54ae4b9e44b928230774ca5620d6d/parameters.p' ... done 00:00:00.49
identifying kmers suitable for use as primers in all 3 ingroup genome sequences
    getting shared ingroup kmers that appear once in each genome ... done 00:01:36.15
    dumping shared kmers to 'primerforge_46f54ae4b9e44b928230774ca5620d6d/sharedKmers.p' ... done 00:00:11.22
    evaluating 2430140 kmers ... done 00:01:51.73
    identified 30413 candidate kmers
    dumping candidate kmers to 'primerforge_46f54ae4b9e44b928230774ca5620d6d/candidates.p' ... done 00:00:01.55
done 00:03:42.94
identifying pairs of primers found in all ingroup sequences ... done 00:00:11.62
    identified 16050 primer pairs shared in all ingroup sequences
    dumping unfiltered pairs to 'primerforge_46f54ae4b9e44b928230774ca5620d6d/pairs.p' ... done 00:00:00.56
removing primer pairs present in the outgroup sequences
    getting outgroup PCR products ... done 00:00:01.03
    filtering primer pairs ... done 00:00:00.54
    processing outgroup results ... done 00:00:00.54
    removed 5905 pairs present in the outgroup (10145 pairs remaining)
    dumping filtered pairs to 'primerforge_46f54ae4b9e44b928230774ca5620d6d/pairs_noOutgroup.p' ... done 00:00:00.54
validating primer pairs with isPcr ... done 00:00:03.44
    removed 3659 pairs not validated by isPcr (6486 pairs remaining)
    dumping validated pairs to 'primerforge_46f54ae4b9e44b928230774ca5620d6d/pairs_noOutgroup_validated.p' ... done 00:00:00.54
writing 6486 primer pairs to 'results.tsv' ... done 00:00:11.62

total runtime: 00:04:14.50
```

As we can see, `primerForge` found 30,413 kmers that were suitable for use as a primer in all three ingroup genomes. It then went on to identify 16,050 primer pairs that would produce PCR products between 500bp and 2000bp in the ingroup genomes. Next, it found that of those 16,050 pairs, 5,905 of them formed PCR products between 300bp and 3000bp in one or more of the outgroup genomes. Finally, it used `isPcr` to validate the remaining 10,145 primer pairs resulting in 6,486 primer pairs being written to file.

### Examining the results
`primerForge` generated `results.tsv`, `primers.bed`, and `primerForge.log`.

#### `results.tsv`
Here are a few lines from `results.tsv`:

|fwd_seq|fwd_Tm|fwd_GC|rev_seq|rev_Tm|rev_GC|i1.gbff_contig|i1.gbff_length|i2.gbff_contig|i2.gbff_length|i3.gbff_contig|i3.gbff_length|o1.gbff_contig|o1.gbff_length|o2.gbff_contig|o2.gbff_length|
|:------|:-----|:-----|:------|:-----|:-----|:------------:|:------------:|:------------:|:------------:|:------------:|:------------:|:------------:|:------------:|:------------:|:------------:|
|AGAAGCAAAGGATATGGGAACAAC|57.1|41.7|AAATCAACACCCTCAATAAGCTCC|57.1|41.7|NZ_LAUX01000078.1|804|NZ_LAUV01000035.1|804|NZ_LAUY01000092.1|804|NA|0|NA|0|
|TCCATCTAATGAAGATCAACCAGG|55.9|41.7|CCCTAATTGTGATGAGTTGACAAC|55.9|41.7|NZ_LAUX01000010.1|710|NZ_LAUV01000042.1|713|NZ_LAUY01000078.1|710|NA|0|NA|0|
|CATCAGCTGTTGTAAATAACCCAC|56.2|41.7|GTGGAGCTATGAAACCAATATCAG|55.3|41.7|NZ_LAUX01000117.1|1694|NZ_LAUV01000064.1|1694|NZ_LAUY01000018.1|1694|NZ_LS483503.1|3212|NA|0|

The first six columns show the forward and reverse sequences (5' --> 3') as well as their melting temperatures and their G+C content (mol %). Next, for each genome it lists the contig and the PCR product size that is predicted be produced by this pair. For example, the first pair of primers are predicted to produce PCR products of 804bp the ingroup genomes, and the binding sites for this primer pair in the files `i1.gbff`, `i2.gbff`, and `i3.gbff` can be found in contigs `NZ_LAUX01000078.1`, `NZ_LAUV01000035.1`, and `NZ_LAUY01000092.1`, respectively. This same pair is not predicted to produce any PCR products in either outgroup genome. Similarly, the third pair is predicted to produce a PCR product size of 1,694bp in all three ingroup genomes, no products `o2.gbff`, and a PCR product >3,000bp in `o1.gbff`.

If a primer pair is predicted to produce multiple products in an outgroup genome, then the contig column and the size column will list contigs and sizes in a comma-separated list linked by position. For example, if a primer pair was expected to produce a product of 1,990bp in `contig_1` and 2,024bp in `contig_2` in the genome file `example.gbff`, then the columns for this genome would look like this:

|example.gbff_contig|example.gbff_length|
|:-----------------:|:-----------------:|
|1990,2024|contig1,contig2|

> [!NOTE]
> Multiple PCR products will only ever be predicted for outgroup genomes as `primerForge` does not allow such primer pairs in the ingroup genome.

#### `primers.bed`
Here are a few lines of `primers.bed`:

```text
NZ_LAUY01000092.1       1420    1444    AGAAGCAAAGGATATGGGAACAAC        0       -
NZ_LAUV01000035.1       11077   11101   AGAAGCAAAGGATATGGGAACAAC        0       +
NZ_LAUX01000078.1       567     591     AGAAGCAAAGGATATGGGAACAAC        0       +
NZ_LAUY01000092.1       640     664     AAATCAACACCCTCAATAAGCTCC        0       +
NZ_LAUV01000035.1       11857   11881   AAATCAACACCCTCAATAAGCTCC        0       -
NZ_LAUX01000078.1       1347    1371    AAATCAACACCCTCAATAAGCTCC        0       -
NZ_LAUY01000057.1       18978   19002   AGTTGGGATTAACCAGACTTCATC        1       +
NZ_LAUV01000024.1       48352   48376   AGTTGGGATTAACCAGACTTCATC        1       +
NZ_LAUX01000061.1       2318    2342    AGTTGGGATTAACCAGACTTCATC        1       +
NZ_LAUY01000057.1       19771   19795   GCTATTTCAAACGCTAAAGCTAGG        1       -
NZ_LAUV01000024.1       49145   49169   GCTATTTCAAACGCTAAAGCTAGG        1       -
NZ_LAUX01000061.1       3111    3135    GCTATTTCAAACGCTAAAGCTAGG        1       -
```

This file is in [BED file format](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bed-format) with one modification: the `score` column (the fifth column) has been overloaded with a "pair number" so that the primer sequences can be linked to one another more easily. For example, pair `0` is `AGAAGCAAAGGATATGGGAACAAC` and `AAATCAACACCCTCAATAAGCTCC` and corresponds with the first entry in `results.tsv`, and pair `1` is `AGTTGGGATTAACCAGACTTCATC` and `GCTATTTCAAACGCTAAAGCTAGG` and corresponds to the second entry in `results.tsv`.


### Finding primer pairs that are target-specific
Let's assume that we want to filter our results to find primer pairs that amplify regions of the _rpoC_ gene in our three ingroup isolates. First, we need to make a file in [BED format](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bed-format) that lists the coordinates of _rpoC_ in our genomes. For _rpoC_ in the ingroup genomes, it would look like this:

```text
NZ_LAUX01000109.1	97	3866	rpoC_i1
NZ_LAUV01000076.1	23	3792	rpoC_i2
NZ_LAUY01000008.1	83	3852	rpoC_i3
```

The first column is the contig, the second column is the start position of the gene, the third column is the end position of the gene (exclusive), and the last column is the name of the entry. We will save this file as `rpoC.bed`.

In order to find the primers that intersect with these genes, we can use `bedtools` to find where our primers intersect with the _rpoC_ gene. If necessary, `bedtools` can be installed with the following command:

```bash
conda install bedtools
```

Once installed, we can use the following command to generate a new bed file that contains only the primers that hit somewhere along the _rpoC_ gene sequences:

```bash
bedtools intersect -a ./primers.bed -b ./rpoC.bed > rpoC_primers.bed
```

The first few lines of the file `rpoC_primers.bed` will look like this:

```text
NZ_LAUY01000008.1       1460    1484    GTTTGAAATATTACCACGAGCTCC        3       +
NZ_LAUV01000076.1       1400    1424    GTTTGAAATATTACCACGAGCTCC        3       +
NZ_LAUX01000109.1       1474    1498    GTTTGAAATATTACCACGAGCTCC        3       +
NZ_LAUY01000008.1       2962    2986    CACCAGACATTCGTCCAATTATTC        3       -
NZ_LAUV01000076.1       2902    2926    CACCAGACATTCGTCCAATTATTC        3       -
NZ_LAUX01000109.1       2976    3000    CACCAGACATTCGTCCAATTATTC        3       -
NZ_LAUY01000008.1       3780    3804    GATCATGAACGAATTGTATCAGGG        74      +
NZ_LAUV01000076.1       3720    3744    GATCATGAACGAATTGTATCAGGG        74      +
NZ_LAUX01000109.1       3794    3818    GATCATGAACGAATTGTATCAGGG        74      +
NZ_LAUY01000008.1       574     598     GGTGTAAATCAATAGCTCCTTCAG        92      +
NZ_LAUV01000076.1       514     538     GGTGTAAATCAATAGCTCCTTCAG        92      +
NZ_LAUX01000109.1       588     612     GGTGTAAATCAATAGCTCCTTCAG        92      +
NZ_LAUY01000008.1       2548    2572    GACTTCAAGATCATGAGTATGCTG        92      -
NZ_LAUV01000076.1       2488    2512    GACTTCAAGATCATGAGTATGCTG        92      -
NZ_LAUX01000109.1       2562    2586    GACTTCAAGATCATGAGTATGCTG        92      -
```

As we can see, the forward and reverse primers for pair numbers 3 and 92 intersect with the _rpoC_ gene in all three isolates. However, only the primer on the (+) strand intersects with _rpoC_ for pair number 74.

## Common error messages (and possible solutions)
### `detected wildcards that are not enclosed in quotes`
This error occurs if you specify a wildcard representing input files without enclosing them in quotes. For example, this will cause the error:

```bash
primerForge --ingroup ./i*gbff
```

and this will fix it:

```bash
primerForge --ingroup "./i*gbff"
```

The same holds true for the `--outgroup` flag.

This error can also occur if you have inadvertently included a space in any of the arguments passed to other flags.

### `invalid or missing file(s)`
This error occurs if the specified file(s) cannot be found or if the file format does not match the `--format` flag (default = `genbank`). Check that the file path is correct and the files can be read. If they are correct, then double check that you have specified `--format fasta` if you are working with `fasta` files.

### `failed to identify a set of kmers shared between the ingroup genomes`
This error occurs if `primerForge` cannot find kmers that are shared in all the ingroup genomes. This can occur if the input genomes are too distantly related, or if one or more of the genomes is of very poor quality. To diagnose this, try repeating the command but include the `--debug` flag. This will report which ingroup genome is causing the number of shared kmers to drop to zero in the `primerForge.log` file.

### `failed to find primer pairs that are absent in the outgroup`
This error occurs if all the primer pairs `primerForge` identified cannot be used to distinguish the ingroup from the outgroup. This most often occurs because one or more members of the outgroup is too closely-related to the ingroup. To diagnose this, try repreating the command but include the `--debug` flag. This will report which outgroup genome is causing the number of shared kmers to drop to zero in the `primerForge.log` file. Alternatively, you can expand your search by widening the ranges passed to the flags `--pcr_prod` and/or `--bad_sizes`.

### `ImportError: /lib64/libstdc++.so.6: version 'GLIBCXX_3.4.20' not found (required by khmer/_khmer.cpython-311-x86_64-linux-gnu.so)`
This error will occur if the C++ library on your machine does not include the symbol `GLIBCXX_3.4.20` which is used by the `khmer` package. This can be resolved by installing the GNU Standard C++ Library using the following `conda` command:

```bash
conda install libstdcxx-ng
```

## Contributing
Thank you for your interest in contributing! Please see the [contributing guidelines](.github/CONTRIBUTING.md) for more information.
