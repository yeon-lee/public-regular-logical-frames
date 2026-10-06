"""Run a contiguous group of one-CPU target shards in one cluster allocation."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--start',type=int,required=True);ap.add_argument('--count',type=int,default=16)
ap.add_argument('--shards',type=int,required=True);ap.add_argument('--seconds',type=float,default=7100);ap.add_argument('--split-weight',type=int,default=4)
ap.add_argument('--node-cap',type=int,default=10**15);ap.add_argument('--output',default='results');args=ap.parse_args()
def run(shard):
    return subprocess.run([sys.executable,str(HERE/'worker.py'),'--shard',str(shard),'--shards',str(args.shards),'--seconds',str(args.seconds),
                           '--split-weight',str(args.split_weight),'--node-cap',str(args.node_cap),'--output',args.output]).returncode
ids=list(range(args.start,min(args.start+args.count,args.shards)))
with ThreadPoolExecutor(max_workers=args.count) as pool:codes=list(pool.map(run,ids))
raise SystemExit(max(codes,default=0))
