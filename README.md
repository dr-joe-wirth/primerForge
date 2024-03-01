# primerForge
software to identify primers that can distinguish genomes

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
                GC-->Tm-->homo-->hair
            end
        end
    end

    %% connections up to candidate kmers
    sharedKmers --> B
    dump1[/"dump to file"/]
    sharedKmers --> dump1
    candidates(["unique, shared kmers; one per start position"])
    hair --> candidates

    %% get primer pairs
    subgraph C["for each genome"]
        bin1["bin overlapping kmers (64bp max)"]
        bin2["get bin pairs"]
        candPair(["candidate primer pairs"])
        sharePair(["shared primer pairs"])

        %% evaluate one kmer pair
        subgraph C0["for each bin pair"]
            size{"is PCR
            size ok?"}
            subgraph C4["for each primer pair"]
                prime{"is 3' end
                G or C?"}
            end
            size --> C4
        end


        %% get shared primer pairs
        subgraph C2["for each candidate primer pair"]
            subgraph C3["for each other genome"]
                pcr{"PCR size ok?"}
            end
        end

        bin1 --> bin2
        bin2 --> C0
        prime --> candPair
        candPair --> C2
        pcr --> sharePair
    end

    allSharePair(["all shared primer pairs"])
    dump2[/"dump to file"/]
    dump3[/"dump to file"/]

    candidates --> dump2
    candidates --> C
    sharePair --> allSharePair
    allSharePair --> dump3

    %% outgroup removal
    outgroup[/"outgroup genomes"/]

    allSharePair --> D0
    outgroup --> D
    subgraph D["for each outgroup genome"]
        subgraph D0["for each primer pair"]
            ogsize{"PCR size outside
            disallowed range?"}
        end
    end
    
    allPairs(["all suitable primer pairs"])
    ogsize --> allPairs

    %% one pair per bin pair
    subgraph E["for each bin pair"]
        keep["keep only one primer
        pair per bin pair"]
    end

    allPairs --> E

    final(["final set of pairs"])
    keep --> final

    write[/"write pairs to file"/]
    plots[/"make plots"/]

    final --> write
    final --> plots
```