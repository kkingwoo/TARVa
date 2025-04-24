import sys
from gatkPipe import GATKPipe
import json
from datetime import datetime
from concurrent import futures

def init_pipe(fis,ref,ds_var,know):
    final_results = []
    info_list = None
    with open(fis, 'r') as fil:
        info_list = json.load(fil)
    with futures.ProcessPoolExecutor(max_workers=12) as mst:
        print('************ Starting parallel processing: FULL GATK PIPELINE for TARVa --> --> ',datetime.now())
        wait_for = [mst.submit(GATKPipe.run_pipeline,inf,ref,ds_var,know) for inf in info_list]
        for fu in futures.as_completed(wait_for):
            current = fu.result()
            final_results.append(current)
    for fi in final_results:
        print(fi)
    return

if __name__=="__main__":
    files = sys.argv[1]
    fast = sys.argv[2]
    dsvar = sys.argv[3]
    known = sys.argv[4]
    run = init_pipe(files,fast,dsvar,known)
