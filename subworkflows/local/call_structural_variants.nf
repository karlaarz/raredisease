//
// A nested subworkflow to call structural variants.
//

include { CALL_SV_MANTA                  } from './call_sv_manta'
include { CALL_SV_MT                     } from './call_sv_MT'
include { CALL_SV_TIDDIT                 } from './call_sv_tiddit'
include { SVDB_MERGE                     } from '../../modules/nf-core/svdb/merge/main'
include { CALL_SV_GERMLINECNVCALLER      } from './call_sv_germlinecnvcaller'
include { CALL_SV_CNVNATOR               } from './call_sv_cnvnator'
include { TABIX_TABIX                    } from '../../modules/nf-core/tabix/tabix/main'

workflow CALL_STRUCTURAL_VARIANTS {

    take:
        ch_genome_bam          // channel: [mandatory] [ val(meta), path(bam) ]
        ch_genome_bai          // channel: [mandatory] [ val(meta), path(bai) ]
        ch_genome_bam_bai      // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_mt_bam_bai          // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_mtshift_bam_bai     // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_bwa_index           // channel: [mandatory] [ val(meta), path(index)]
        ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai          // channel: [mandatory] [ val(meta), path(fai) ]
        ch_mtshift_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_case_info           // channel: [mandatory] [ val(case_info) ]
        ch_target_bed          // channel: [mandatory for WES] [ val(meta), path(bed), path(tbi) ]
        ch_genome_dictionary   // channel: [optional; used by mandatory for GATK's cnvcaller][ val(meta), path(dict) ]
        ch_svcaller_priority   // channel: [mandatory] [ val(["var caller tag 1", ...]) ]
        ch_readcount_intervals // channel: [optional; used by mandatory for GATK's cnvcaller][ path(intervals) ]
        ch_ploidy_model        // channel: [optional; used by mandatory for GATK's cnvcaller][ path(ploidy_model) ]
        ch_gcnvcaller_model    // channel: [optional; used by mandatory for GATK's cnvcaller][ path(gcnvcaller_model) ]

    main:
        ch_versions = Channel.empty()
        ch_merged_svs = Channel.empty()
        ch_merged_tbi = Channel.empty()

                if (!params.analysis_type.equals("mito")) {
            CALL_SV_MANTA (ch_genome_bam, ch_genome_bai, ch_genome_fasta, ch_genome_fai, ch_case_info, ch_target_bed)
            ch_versions = ch_versions.mix(CALL_SV_MANTA.out.versions)

            manta_vcf = CALL_SV_MANTA.out.filtered_diploid_sv_vcf_tbi
                .map { meta, vcf -> [ meta.id, vcf ] }
                .groupTuple()
        }

        if (params.analysis_type.equals("wgs")) {
            CALL_SV_TIDDIT (ch_genome_bam_bai, ch_genome_fasta, ch_bwa_index, ch_case_info)
            ch_versions = ch_versions.mix(CALL_SV_TIDDIT.out.versions)

            tiddit_vcf = CALL_SV_TIDDIT.out.vcf
                .map { meta, vcf -> [ meta.id, vcf ] }
                .groupTuple()

            CALL_SV_CNVNATOR (ch_genome_bam_bai, ch_genome_fasta, ch_genome_fai, ch_case_info)
            ch_versions = ch_versions.mix(CALL_SV_CNVNATOR.out.versions)

            cnvnator_vcf = CALL_SV_CNVNATOR.out.vcf
                .map { meta, vcf -> [ meta.id, vcf ] }
                .groupTuple()
        }

        if (!(params.skip_tools && params.skip_tools.split(',').contains('germlinecnvcaller'))) {
            CALL_SV_GERMLINECNVCALLER (ch_genome_bam_bai, ch_genome_fasta, ch_genome_fai, ch_readcount_intervals, ch_genome_dictionary, ch_ploidy_model, ch_gcnvcaller_model, ch_case_info)
            ch_versions = ch_versions.mix(CALL_SV_GERMLINECNVCALLER.out.versions)

            gcnvcaller_vcf = CALL_SV_GERMLINECNVCALLER.out.genotyped_filtered_segments_vcf
                .map { meta, vcf -> [ meta.id, vcf ] }
                .groupTuple()
        }

        // merge — agrupado por case_id
        if (params.skip_tools && params.skip_tools.split(',').contains('germlinecnvcaller')) {
            if (params.analysis_type.equals("wgs")) {
                tiddit_vcf
                    .join(manta_vcf)
                    .join(cnvnator_vcf)
                    .map { case_id, tiddit, manta, cnvnator -> [ case_id, tiddit + manta + cnvnator ] }
                    .set { vcf_list }
            } else if (!params.analysis_type.equals("mito")) {
                manta_vcf
                    .set { vcf_list }
            }
        } else if (params.analysis_type.equals("wgs")) {
            tiddit_vcf
                .join(manta_vcf)
                .join(gcnvcaller_vcf)
                .join(cnvnator_vcf)
                .map { case_id, tiddit, manta, gcnv, cnvnator -> [ case_id, tiddit + manta + gcnv + cnvnator ] }
                .set { vcf_list }
        } else if (!params.analysis_type.equals("mito")) {
            manta_vcf
                .join(gcnvcaller_vcf)
                .map { case_id, manta, gcnv -> [ case_id, manta + gcnv ] }
                .set { vcf_list }
        }

        if (!params.analysis_type.equals("mito")) {
            ch_case_info
                .map { case_info -> [ case_info.id, case_info ] }
                .join(vcf_list)
                .map { case_id, case_info, vcfs -> [ case_info, vcfs ] }
                .set { merge_input_vcfs }

            SVDB_MERGE (merge_input_vcfs, ch_svcaller_priority, true)
            TABIX_TABIX (SVDB_MERGE.out.vcf)
            ch_merged_svs = SVDB_MERGE.out.vcf
            ch_merged_tbi = TABIX_TABIX.out.tbi
            ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
            ch_versions = ch_versions.mix(SVDB_MERGE.out.versions)
        }


    emit:
        vcf      = ch_merged_svs // channel: [ val(meta), path(vcf)]
        tbi      = ch_merged_tbi // channel: [ val(meta), path(tbi)]
        versions = ch_versions   // channel: [ path(versions.yml) ]
}
