import sys
import os
import re
import pickle
import argparse
import pandas as pd
from natsort import natsorted
from statistics import median


class Logger(object):
    def __init__(self, save_dir, source):
        self.terminal = sys.stdout
        self.log = open(f'{save_dir}/{source}_norm_log.log', 'a')
   
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        pass   


def run (CSV, mfname, save_dir):
    print("Starting normalization")
    
    source_df = pd.read_csv(CSV)
    source = os.path.basename(os.path.normpath(CSV))
    source = re.sub(r'(.*)\.[cC][sS][vV]$', r'\1', source)
    sys.stderr = Logger(save_dir, source) 
    sys.stdout = Logger(save_dir, source) 

    me_markers = pd.read_csv(mfname)
    me_markers = {col: me_markers[col].dropna().to_list() for col in me_markers.columns}
    
    for scene in natsorted(set(source_df.scene)):
        print (source, scene)
        tmp_scene = source_df[source_df.scene==scene]
    #     tmp_scene = tmp_scene.drop(tmp_scene.iloc[:,0:7], axis=1)
        tmp_scene = tmp_scene.drop(['scene'], axis=1)
        Min = tmp_scene.min()
 
        for pos in list(me_markers):
            try:
               nf = pickle.load(open(f'{save_dir}/thresh_dicts/{source}/{source}_{pos}_thresh_dict.pkl', "rb"))
            except:
               print('error opening', f'{save_dir}/thresh_dicts/{source}/{source}_{pos}_thresh_dict.pkl')
               del me_markers[pos]
               continue
            nf = [nf[source][scene][pos]]
            nf = [ele for ele in ({key: val for key, val in sub.items() if val} for sub in nf) if ele]
            nf = list(nf[0].values())
            nf = [element for sublist in nf for element in sublist]
            nf = median(nf)
            source_df.loc[source_df.scene == scene,pos] = (source_df.loc[source_df.scene == scene,pos] - Min[pos]) / (nf - Min[pos])
    source_df = pd.concat([source_df['scene'], source_df[me_markers.keys()]], axis=1)
    source_df.to_csv(f'{save_dir}/{source}_RESTORE.csv', index=False)
        
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='RESTPRE normalization')
    parser.add_argument('--CSV', type=str, default=None, metavar='S', help='which CSV file to process')
    parser.add_argument('--mfname', type=str, default=None, metavar='S', help='which marker file to process')
    parser.add_argument('--save_dir', type=str, default='~/', metavar='S', help='same dir as calculate_threshols.py')
    args = parser.parse_args()
    
    run(args.CSV, args.mfname, args.save_dir)

