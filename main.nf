#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Automated pipeline for Freyja analysis of SRA data
 */

input = file(params.input)

ref = file(params.ref)
bedfiles = file(params.bedfiles)
baseDir = file("$baseDir")
annot = file(params.annot)

// Import modules
include {
    GET_ACCESSIONS;
    GET_AMPLICON_SCHEME;
    SRA_PREFETCH;
    FASTERQ_DUMP;
} from "./modules/sra.nf"

include {
    CUTADAPT_TRIM;
    MINIMAP2;
    MINIMAP2_UNKNOWN_PRIMER;
    SAMTOOLS_1;
    SAMTOOLS_2;
    IVAR_TRIM;
} from "./modules/preprocessing.nf"

include {
    FREYJA_VARIANTS;
    FREYJA_DEMIX;
    FREYJA_COVARIANTS;
} from "./modules/freyja.nf"

include {
    AGGREGATE_VARIANTS
    AGGREGATE_DEMIX;
    AGGREGATE_COVARIANTS;
} from "./modules/aggregate.nf"

workflow fetch_sra {

    Channel
        .fromPath(input)
        .set { input_ch }

    GET_ACCESSIONS(input_ch)
        .splitCsv()
        .map { line -> line.join('') }
        .take(params.num_samples)
        .set { acc_ch }


    GET_AMPLICON_SCHEME(acc_ch, input)
        .set { primer_scheme_ch }

    SRA_PREFETCH(primer_scheme_ch)

    FASTERQ_DUMP(SRA_PREFETCH.out)
        .branch {
            unknown_primer: it[1].text == 'unknown'
            known_primer: it[1].text != 'unknown'
        }
        .set { fq_ch }

    process_unknown_primer(fq_ch.unknown_primer)
    process_known_primer(fq_ch.known_primer)

}

workflow process_unknown_primer {
    take: unknown_primer_fastq_ch

    main:
    CUTADAPT_TRIM(unknown_primer_fastq_ch)
    MINIMAP2_UNKNOWN_PRIMER(CUTADAPT_TRIM.out, ref)
    SAMTOOLS_1(MINIMAP2_UNKNOWN_PRIMER.out)

    freyja(SAMTOOLS_1.out)
}

workflow process_known_primer {
    take: known_primer_fastq_ch

    main:
    known_primer_fastq_ch
        .map { it[1].text }
        .set { primer_ch }
    MINIMAP2(known_primer_fastq_ch, ref)
    SAMTOOLS_1(MINIMAP2.out)
    IVAR_TRIM(SAMTOOLS_1.out, primer_ch, bedfiles)
    SAMTOOLS_2(IVAR_TRIM.out)

    freyja(SAMTOOLS_2.out)
}

workflow freyja {
    take:
    sra_accession
    input_bam
    bam_index
    
    main:
    FREYJA_VARIANTS(sra_accession, input_bam, bam_index, ref)
        .collect()
        .map { it[1] }
        .set { variants_ch }

    FREYJA_DEMIX(FREYJA_VARIANTS.out, params.eps, params.depthCutoff)
        .collect()
        .set { demix_ch }

    // FREYJA_COVARIANTS(sra_accession, input_bam, bam_index, ref, annot)
    //     .collect()
    //     .set { covariants_ch }
}


workflow rerun_demix {
    Channel
        .fromFilePairs("${params.input_dir}/SRR*{variants,depths}.tsv")
        .set { variants_ch }

    FREYJA_DEMIX(variants_ch)
        .collect()
        .set { demix_ch }
}
