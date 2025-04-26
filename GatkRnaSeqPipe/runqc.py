import os
import subprocess

class RunQC:
    @staticmethod
    def mono_50_trim(i1,i2,op1,oup1,op2,oup2,adf):
        t_comm = f"trimmomatic PE -threads 4 {i1} {i2} {op1} {oup1} {op2} {oup2} ILLUMINACLIP:{adf}:2:30:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:35"
        try:
            subprocess.run(t_comm, shell=True, check=True)
            print(f"Trimmomatic completed successfully for {i1} and {i2}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Trimmomatic failed with error: {e}")
            return False

    @staticmethod
    def mono_76_trim(i1,i2,op1,oup1,op2,oup2,adf):
        t_comm = f"trimmomatic PE -threads 4 {i1} {i2} {op1} {oup1} {op2} {oup2} ILLUMINACLIP:{adf}:2:30:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:40"
        try:
            subprocess.run(t_comm, shell=True, check=True)
            print(f"Trimmomatic completed successfully for {i1} and {i2}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Trimmomatic failed with error: {e}")
            return False

    @staticmethod
    def bulk_150_trim(i1,i2,op1,oup1,op2,oup2,adf):
        t_comm = f"trimmomatic PE -threads 8 {i1} {i2} {op1} {oup1} {op2} {oup2} ILLUMINACLIP:{adf}:2:30:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:50"
        try:
            subprocess.run(t_comm, shell=True, check=True)
            print(f"Trimmomatic completed successfully for {i1} and {i2}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Trimmomatic failed with error: {e}")
            return False

    @staticmethod
    def run_qc(diction,length,tissue,adapter_file):
        trimmed_50,trimmed_76,trimmed_150 = [50],[76],[150]
        for key in diction.keys():
            for ri in diction[key].keys():
                files = diction[key][ri]
                if len(files) == 2:
                    f1,f2 = files[0],files[1]
                    pa1,pa2 = f1.split('/'),f2.split('/')
                    pat1,pat2 = '/'.join(pa1[:-1])+'/Trimmed/','/'.join(pa2[:-1])+'/Trimmed/'
                    if tissue == "monocytes":
                        pf1s,upf1s,pf2s,upf2s = ri+'_paired_'+pa1[-1].split('_')[2],ri+'_unpaired_'+pa1[-1].split('_')[2],ri+'_paired_'+pa2[-1].split('_')[2],ri+'_unpaired_'+pa2[-1].split('_')[2]
                        pap1,paup1,pap2,paup2 = pat1+pf1s,pat1+upf1s,pat2+pf2s,pat2+upf2s
                        if length == '50':
                            trimmed_50.append((pap1,pap2))
                            #trim_50 = RunQC.mono_50_trim(f1,f2,pap1,paup1,pap2,paup2,adapter_file)
                            #trimmed_50.append(trim_50)
                        elif length == '76':
                            trimmed_76.append((pap1,pap2))
                            #trim_76 = RunQC.mono_76_trim(f1,f2,pap1,paup1,pap2,paup2,adapter_file)
                            #trimmed_76.append(trim_76)
                        else:
                            trimmed_76.append(f"{length} was a length not factored in monocyte dataset")
                    if tissue == "bulkbrain":
                        check1,check2 = pa1[-1].split('_'), pa2[-1].split('_')   
                        pf1s,upf1s,pf2s,upf2s = ri+'_paired_'+'_'.join(check1[-2:]),ri+'_unpaired_'+'_'.join(check1[-2:]),ri+'_paired_'+'_'.join(check2[-2:]), ri+'_unpaired_'+'_'.join(check2[-2:])
                        pap1,paup1,pap2,paup2 = pat1+pf1s,pat1+upf1s,pat2+pf2s,pat2+upf2s
                        if length == '150':
                            trimmed_150.append((pap1,pap2))
                            #trim_150 = RunQC.bulk_150_trim(f1,f2,pap1,paup1,pap2,paup2,adapter_file)
                            #trimmed_150.append(trim_150)
                        else:
                            trimmed_150.append("{length} was a length not factored in bulk_brain dataset")
                            
                    
        return trimmed_50,trimmed_76,trimmed_150


