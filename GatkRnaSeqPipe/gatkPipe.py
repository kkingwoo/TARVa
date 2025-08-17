import os
import subprocess
import gzip

class GATKPipe:
    @staticmethod
    def check_bam(b):
        if os.path.exists(b):
            try:
                checks = f"samtools quickcheck -vvv {b}"
                subprocess.run(checks,check=True,shell=True)
                return True
            except subprocess.CalledProcessError:
            # If samtools quickcheck fails, delete the invalid BAM file
                print(f"Invalid or incomplete BAM detected: {b}. Removing and reprocessing.")
                os.remove(b)
            except Exception as e:
                print(f"Unexpected error during validation: {e}")
        return False

    @staticmethod
    def fastq2sam(f1,f2,rid,ds):
        name = None
        checks = None
        with gzip.open(f1,'rt') as fq:
            name = next(line.strip() for line in fq if line.startswith('@'))
        n = name.split(':')[1]
        na = name.split(":")[0][1:]
        rg = f"{na}.{n}"
        pu = rg.split(".")[-1]
        lb = name.split(':')[-1]
        d = os.path.dirname(f1)
        obam = f"{d}/{rid}_varcall_unaln.bam"
        if GATKPipe.check_bam(obam):
            return f"BAM is valid and complete: {obam}"
        if not os.path.exists(obam):
        comm1 = f"gatk FastqToSam -F1 {f1} -F2 {f2} -O {obam} -SM {rid} -RG {rg} -DS {ds}"
        subprocess.run(comm1,check=True,shell=True)
        return obam,d,rg,pu,lb
    
    @staticmethod
    def samsort1(ibam,di,rid):
        out_srtd = f"{di}/{rid}_varcall_srtByQN.bam"
        check_state = f"Everything is fine for BAM: {out_srtd}"
        comm2 = f'gatk SortSam -I {ibam} -O {out_srtd} -SO queryname'
        if GATKPipe.check_bam(out_srtd):
            return f"BAM is valid and complete: {out_srtd}"
        if not os.path.exists(out_srtd):
        subprocess.run(comm2,check=True,shell=True)
        return out_srtd

    @staticmethod
    def samsort2(ibam,di,rid):
        mrk_dupes = f"{di}/{rid}_marked_dupes_srt.bam"
        if GATKPipe.check_bam(mrk_dupes):
            return f"BAM is valid and complete: {mrk_dupes}"
        if not os.path.exists(mrk_dupes):
        comm6 = f"gatk SortSam -I {ibam} -O {mrk_dupes} --SORT_ORDER coordinate"
        subprocess.run(comm6,check=True,shell=True)
        return mrk_dupes

    @staticmethod
    def mergebamalignment(unsrt_bm,srt_bm,dirs,rid,ref):
        merged_bam = f"{dirs}/{rid}_varcall_merge.bam"
        if GATKPipe.check_bam(merged_bam):
            return f"BAM is valid and complete: {merged_bam}"
        if not os.path.exists(merged_bam):
        comm3 = f"gatk MergeBamAlignment -ALIGNED {srt_bm} -UNMAPPED {unsrt_bm} -O {merged_bam} -R {ref}"
        subprocess.run(comm3,check=True,shell=True)
        return merged_bam

    @staticmethod
    def samtools(mrgd,di,rid):
        srtd = f"{di}/{rid}_srt.bam"
        if GATKPipe.check_bam(srtd):
            return f"BAM is valid and complete: {srtd}"
        if not os.path.exists(srtd):
        comm4 = f"samtools sort {mrgd} -o {srtd}"
        subprocess.run(comm4,check=True,shell=True)
        return srtd

    @staticmethod
    def markdupes(srt,ref,di,rid):
        marked = f"{di}/{rid}_markedDupes.bam"
        metrics = f"{di}/{rid}_metrics.txt"
        if GATKPipe.check_bam(marked):
            return f"BAM is valid and complete: {marked}"
        if not os.path.exists(marked):
            if os.path.exists(metrics):
                os.remove(metrics)
        comm5 = f"gatk MarkDuplicates -I {srt} -O {marked} -M {metrics} --TAGGING_POLICY All -R {ref}"
        subprocess.run(comm5,check=True,shell=True)
        return marked 

    @staticmethod
    def splitncigar(inbam,di,rid,ref):
        split = f"{di}/{rid}_split_Ncigar_srt.bam"
        if GATKPipe.check_bam(split):
            return f"BAM is valid and complete: {split}"
        if not os.path.exists(split):
        comm7 = f"gatk SplitNCigarReads -R {ref} -I {inbam} -O {split}"
        subprocess.run(comm7,check=True,shell=True)
        return split

    @staticmethod
    def addreplacereadgroups(inbam,di,rid,rgs,dss,lbs,pus):
        add_replace = f"{di}/{rid}_add_replace.bam"
        if GATKPipe.check_bam(add_replace):
            return f"BAM is valid and complete: {add_replace}"
        if not os.path.exists(add_replace):
        comm8 = f"gatk AddOrReplaceReadGroups -I {inbam} -O {add_replace} --RGSM {rid} --RGID {rgs} --RGDS {dss} --RGPL illumina --RGLB {lbs} --RGPU {pus}"
        subprocess.run(comm8,check=True,shell=True)
        return add_replace

    @staticmethod
    def baserecal(inbam,di,rid,ref,known):
        base_recal = f"{di}/{rid}_base_recal.table"
        comm9 = f"gatk BaseRecalibrator -I {inbam} -R {ref} --known-sites {known} -O {base_recal}"
        subprocess.run(comm9,check=True,shell=True)
        return base_recal

    @staticmethod
    def applybqsr(inbam,intab,di,rid,ref):
        bqsr_bam = f"{di}/{rid}_bqsr.bam"
        if GATKPipe.check_bam(bqsr_bam):
            return f"BAM is valid and complete: {bqsr_bam}"
        if not os.path.exists(bqsr_bam):
        comm10 = f"gatk ApplyBQSR -R {ref} -I {inbam} --bqsr-recal-file {intab} -O {bqsr_bam}"
        subprocess.run(comm10,check=True,shell=True)
        return bqsr_bam

    @staticmethod
    def analyzecovariates(intable,di,rid):
        out_pdf = f"{di}/{rid}_covars_analyzed.pdf"
        comm11 = f"gatk AnalyzeCovariates -bqsr {intable} -plots {out_pdf}"
        subprocess.run(comm11,check=True,shell=True)
        return out_pdf

    @staticmethod
    def bam_roundup(b1,b2,b3):
        for b in [b1,b2,b3]:
            if os.path.exists(b):
                comm12 = f"rm {b}"
                subprocess.run(comm12,check=True,shell=True)
        return

    @staticmethod
    def haplocall(inbam,di,rid,ref):
        out_vcf = f"{di}/{rid}_haplo.vcf"
        comm13 = f"gatk HaplotypeCaller -I {inbam} -O {out_vcf} -R {ref} --output-mode EMIT_ALL_CONFIDENT_SITES"
        if not os.path.exists(out_vcf):
        subprocess.run(comm13,check=True,shell=True)
        return out_vcf

    @staticmethod
    def filter_vars(invcf,di,rid,ref):
        filt_vcf = f"{di}/{rid}_filtered.vcf"
        comm14 = f"gatk VariantFiltration -R {ref} -V {invcf} -O {filt_vcf} --filter-name qual_filt --filter-expression \"QUAL < 25.0\" --filter-name qual_depth_filt --filter-expression \"QD < 10.0\" --filter-name depth_filt --filter-expression \"DP < 20\""
        subprocess.run(comm14,check=True,shell=True)
        return filt_vcf

    @staticmethod
    def run_pipeline(fi_set,ref_fa,ds,known):
        f1,f2,bam,ri = fi_set[0],fi_set[1],fi_set[2],fi_set[3]
        unalign_bam = GATKPipe.fastq2sam(f1,f2,ri,ds)
        unaln_bam,dirs,rg,pu,lb = unalign_bam[0],unalign_bam[1],unalign_bam[2],unalign_bam[3],unalign_bam[4]
        sorted_bam = GATKPipe.samsort1(bam,dirs,ri)
        merged_bam = GATKPipe.mergebamalignment(unaln_bam,sorted_bam,dirs,ri,ref_fa)
        sorted_out = GATKPipe.samtools(merged_bam,dirs,ri)
        mrkd_out = GATKPipe.markdupes(sorted_out,ref_fa,dirs,ri)
        mrkd_srt_out = GATKPipe.samsort2(mrkd_out,dirs,ri)
        splitn_out = GATKPipe.splitncigar(mrkd_srt_out,dirs,ri,ref_fa)
        add_replace_rg = GATKPipe.addreplacereadgroups(splitn_out,dirs,ri,rg,ds,lb,pu)
        baserec = GATKPipe.baserecal(add_replace_rg,dirs,ri,ref_fa,known)
        apply_bsqr = GATKPipe.applybqsr(add_replace_rg,baserec,dirs,ri,ref_fa)
        analyze_covars = GATKPipe.analyzecovariates(baserec,dirs,ri)
        bams_removed = GATKPipe.bam_roundup(unaln_bam,merged_bam,mrkd_out)
        haplos = GATKPipe.haplocall(apply_bsqr,dirs,ri,ref_fa)
        #haplos = f"{dirs}/{ri}_haplo.vcf"
        filtered = GATKPipe.filter_vars(haplos,dirs,ri,ref_fa)
        string = f"Data preprocessing complete; final output file:\n{filtered}"
        return string
