<img src="assets/logo.png" alt="Logo" width="250" height="250">

# primerForge

software to identify primers that can be used to distinguish genomes

## Installation
`primerForge` is incompatible with versions python3.8 and below and python3.12 and above.

### `pip` installation
`primerForge` can be installed with the following commands:

```shell
pip install primerforge
conda install ispcr
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

### Docker Installation

A Docker image for the latest release is available at [DockerHub](https://hub.docker.com/r/jwirth/primerforge)

### Checking installation

If `primerForge` is installed correctly, then the following command should execute without errors:

```shell
primerForge --check_install
```

If you installed manually, you may need to use the following command instead

```shell
python primerForge.py --check_install
```

### Running unit tests
In order to run unit tests, install `primerForge` using [the instructions above](#installation). You will also need to clone the repository if you haven't already:

```bash
git clone https://github.com/dr-joe-wirth/primerForge.git
```

Once installed and cloned, run the following commands to run the unit tests:

> [!NOTE]
> Running `results_test.py` may take up to three hours to complete

```bash
python3 -m unittest discover -s .primerForge/bin/unit_tests/ -p "*_test.py"
```


## Usage

```text
usage:
    primerForge [-ioaubfpgtrdnkvh]

required arguments:
    -i, --ingroup        [file] ingroup filename or a file pattern inside double-quotes (eg."*.gbff")

optional arguments: 
    -o, --out            [file] output filename for primer pair data (default: results.tsv)
    -a, --analysis       [file] output basename for primer analysis data (default: distribution)
    -u, --outgroup       [file(s)] outgroup filename or a file pattern inside double-quotes (eg."*.gbff")
    -b, --bad_sizes      [int,int] a range of PCR product lengths that the outgroup cannot produce (default: same as '--pcr_prod')
    -f, --format         [str] file format of the ingroup and outgroup genbank|fasta (default: genbank)
    -p, --primer_len     [int(s)] a single primer length or a range specified as 'min,max' (default: 16,20)
    -g, --gc_range       [float,float] a min and max percent GC specified as a comma separated list (default: 40.0,60.0)
    -t, --tm_range       [float,float] a min and max melting temp (Tm) specified as a comma separated list (default: 55.0,68.0)
    -r, --pcr_prod       [int(s)] a single PCR product length or a range specified as 'min,max' (default: 120,2400)
    -d, --tm_diff        [float] the maximum allowable Tm difference between a pair of primers (default: 5.0)
    -n, --num_threads    [int] the number of threads for parallel processing (default: 1)
    -k, --keep           keep intermediate files (default: False)
    -v, --version        print the version
    -h, --help           print this message
    --check_install      check installation
    --debug              run in debug mode (default: False)
```

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
    dump1[/"dump to file"/]
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
    prime --> temp --> hetero --> candPair
    candPair --> D
    pcr --> allSharePair

    %% one pair per bin pair
    subgraph E["for each bin pair"]
        keep["keep only one primer pair"]
    end
    
    selectedSharePair(["selected shared primer pairs"])
    dump2[/"dump to file"/]
    dump3[/"dump to file"/]

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
    dump4[/"dump to file"/]

    ispcr["filter pairs using isPcr"]
    final(["final primer pairs"])
    dump5[/"dump to file"/]
    write[/"write pairs to file"/]


    noout --> ispcr
    ispcr --> final
    final --> dump5
    final --> write
```

## Example usage with a test case
### Motivation
In this example, we will use `primerForge` to find primer pairs of size between 18bp and 24bp that can be used to differentiate three strains of _Mycoplasma mycoides_ subspecies mycoides (the "ingroup") from two strains of _Mycoplasma mycoides_ subspecies capri (the "outgroup"). The primer pairs identified by `primerForge` are predicted to produce a single PCR product between 500bp and 2000bp in the ingroup. These same primer pairs are predicted to produce PCR products <300bp, >3000bp, or no PCR products at all.

### preparing the workspace
In order to get started, create a directory called `mycoplasma_test` and download the following _Mycoplasma mycoides_ genomes using the following commands:

```bash
mkdir ./mycoplasma_test
wget -q -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/034/305/GCF_003034305.1_ASM303430v1/GCF_003034305.1_ASM303430v1_genomic.gbff.gz | gzip -d > ./mycoplasma_test/i1.gbff
wget -q -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/034/275/GCF_003034275.1_ASM303427v1/GCF_003034275.1_ASM303427v1_genomic.gbff.gz | gzip -d > ./mycoplasma_test/i2.gbff
wget -q -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/034/345/GCF_003034345.1_ASM303434v1/GCF_003034345.1_ASM303434v1_genomic.gbff.gz | gzip -d > ./mycoplasma_test/i3.gbff
wget -q -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/489/555/GCF_900489555.1_MMC68/GCF_900489555.1_MMC68_genomic.gbff.gz | gzip -d > mycoplasma_test/o1.gbff
wget -q -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/389/745/GCF_018389745.1_ASM1838974v1/GCF_018389745.1_ASM1838974v1_genomic.gbff.gz | gzip -d > ./mycoplasma_test/o2.gbff
```

If you cannot download the genbank files using `wget`, you can download them manually from NCBI by replacing `ftp://` with `http://` and copying and pasting the addresses into your web browser.

Install `primerForge` and confirm it is installed properly [as described above](#installation).

### Running `primerForge`
Run `primerForge` using the following command:

```bash
cd ./mycoplasma_test
primerForge --ingroup "./i*gbff" --outgroup "./o*gbff" --pcr_prod 500,2000 --bad_sizes 300,3000 --primer_len 18,24
```

  * The `--ingroup` and `--outgroup` flags are both file patterns for the ingroup and outgroup genomes, respectively. It is important that this pattern is enclosed in double-quotes as shown above.
  * The `pcr_prod` flag indicates that we want ingroup products between 500 and 2000bp.
  * The `--bad_sizes` flag indicates that we do not want outgroup products between 300bp and 3000bp.
  * The `--primer_len` flag indicates that our primers are allowed to be between 18bp an 24bp.

You can get a list of all available flags using the command `primerForge --help`.

After running that command, you should see something like this printed to the screen:

```text

```

As we can see, `primerForge` found 5615 kmers that were suitable for use as a primer in all three ingroup genomes. It then went on to identify **___** primer pairs that would produce PCR products between 500bp and 2000bp in the ingroup genomes. Next, it found that of those **___** pairs, **___** of them formed PCR products between 300bp and 3000bp in one or more outgroup genomes. Finally, it used `isPcr` to validate the remaining **___** primer pairs resulting in **___** primer pairs being written to file.

`primerForge` generated `results.tsv`, the file that contains the sequences and details for each primer pair, and `primerForge.log`. Here are the first few lines from this file:

```text

```

The first six columns show the forward and reverse sequences (5' --> 3') as well as their melting temperatures and their G+C %. Next, for each genome it lists the contig and the PCR product size(s) that is predicted be produced by this pair.
