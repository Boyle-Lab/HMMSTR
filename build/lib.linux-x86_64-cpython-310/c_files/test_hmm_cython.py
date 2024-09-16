import hmm as hmm
import argparse

#read in inputs for run_main
parser = argparse.ArgumentParser()
parser.add_argument("model_file", type=str)
parser.add_argument("seq_len", type=str)
parser.add_argument("seq", type=str)
parser.add_argument("repeat_len", type=str)
parser.add_argument("prefix_idx", type=str)
args = parser.parse_args()
hmm.run_main(args.model_file,args.seq_len,args.seq,args.repeat_len,args.prefix_idx)

