
import sys
import os
import re
import pickle
import argparse
import warnings
import numpy as np
import pandas as pd
import holoviews as hv
from natsort import natsorted
from data import nested_dict, process_data
from clustering import get_ssc_thresh


warnings.simplefilter('ignore')
hv.extension('matplotlib','bokeh')

np.random.seed(seed=250)

class Logger(object):
    def __init__(self, save_dir, source):
        self.terminal = sys.stdout
        self.log = open(f'{save_dir}/{source}_thresh_log.log', 'a')
   
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        pass   


def run(CSV, mfname, pos_idx, save_dir, floor):
    
    
    source_df = pd.read_csv(CSV)
    source = os.path.basename(os.path.normpath(CSV))
    source = re.sub(r'(.*)\.[cC][sS][vV]$', r'\1', source)
    
    source_threshs_dir = f'{save_dir}/img/threshs/{source}'
    os.makedirs(source_threshs_dir, exist_ok=True)
    sys.stderr = Logger(save_dir, source) 
    sys.stdout = Logger(save_dir, source) 

    me_markers = pd.read_csv(mfname)
    me_markers = {col: me_markers[col].dropna().to_list() for col in me_markers.columns}
    pos = list(me_markers.keys())[int(pos_idx)]

    thresh_dict = nested_dict()

    for neg in me_markers[pos]:
        
        print(pos, neg)

        figs = []

        for scene in natsorted(set(source_df.scene)):
            
            print(source, scene)

            tmp_df = source_df[source_df.scene==scene][[pos,neg]]
            tmp_df = process_data(tmp_df,pos,neg,floor)
            thresh, clusters = get_ssc_thresh(tmp_df)

            thresh_dict[source][scene][pos][neg] = [thresh]

            xlim = (source_df[pos].quantile(0.001),
                    source_df[pos].quantile(0.999))

            ylim = (source_df[neg].quantile(0.001),
                    source_df[neg].quantile(0.999))

            scatters = hv.Overlay([hv.Scatter(i).opts(xlabel=pos,
                                                      ylabel=neg,
                                                      xlim=xlim,
                                                      ylim=ylim,
                                                      alpha=0.1) for i in clusters])

            thresh_line = hv.VLine(thresh).opts(color='black')

            fig = hv.Overlay([scatters, thresh_line],label=f'scene {str(scene)}')

            figs.append(fig)
        
        layout = hv.Layout(figs).opts(title=source,sublabel_format='').cols(4)
        hv.save(layout,f'{source_threshs_dir}/{source}_{pos}_{neg}.png')
    
    thresh_dir = f'{save_dir}/thresh_dicts/{source}'
    os.makedirs(thresh_dir,exist_ok=True)
    
    pickle.dump(thresh_dict, open(f'{thresh_dir}/{source}_{pos}_thresh_dict.pkl','wb'))
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='RESTPRE calculate threshold')
    parser.add_argument('--CSV', type=str, default=None, metavar='S', help='which CSV file to process')
    parser.add_argument('--mfname', type=str, default=None, metavar='S', help='which marker file to process')
    parser.add_argument('--pos_idx', type=str, default=None, metavar='S', help='which marker to normalize')
    parser.add_argument('--save_dir', type=str, default='~/', metavar='S', help='where to save')
    parser.add_argument('--floor', type=float, default=0.0, metavar='S', help='set floor of intensity')
    
    args = parser.parse_args()
    
    run(args.CSV, args.mfname, args.pos_idx, args.save_dir, args.floor)



