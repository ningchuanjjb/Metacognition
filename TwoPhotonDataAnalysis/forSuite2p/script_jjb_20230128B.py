import numpy as np
import sys
# option to import from github folder
#sys.path.insert(0, 'C:/Users/JJB/Anaconda3/envs/suite2p/Lib/site-packages/suite2p')
sys.path.append('C:/Users/JJB/Anaconda3/envs/suite2p/Lib/site-packages/suite2p')
from suite2p import run_s2p
sys.path.append('C:/ASDROOT/STUDY/TwoPhotonDataAnalysis/ForSuite2p')
#from ForSuite2p import mydefault_ops
import mydefault_ops
import glob
import os



current_rootPath = 'C:/ASDROOT/STUDY/twoPhotonData_motionCorrected/113Recording_20230123A_Ding_Site09B_sameFOV0122'


target_dir_prefix = 'Result'


file_list = glob.glob((current_rootPath+'/'+target_dir_prefix+'*'))
print('Find %d dir!\n' %(len(file_list)),end="")
for tempi in range(len(file_list)):
   temp_fullDir = os.path.split(file_list[tempi])
   temp_dir = temp_fullDir[1]
   print('%s\n' %(temp_dir),end="")
   
target_fullDir = os.path.split(file_list[len(file_list)-1])
target_dir = target_fullDir[1]

#target_dir = 'Result20230128T163621'

print('Choose %s to analyse!\n' %(target_dir),end="")
   
   



current_path = current_rootPath + '/' + target_dir

# set your options for running
ops = mydefault_ops.default_ops() # populates ops with the default options

ops['threshold_scaling'] = 0.70 # 1.5-->0.75

#ops['ops.reg_file'] = current_path + '/plane0/data1.bin'

# provide an h5 path in 'h5py' or a tiff path in 'data_path'
# db overwrites any ops (allows for experiment specific settings)
db = {
      'h5py': [], # a single h5 file path
      'h5py_key': 'mov',
      'look_one_level_down': False, # whether to look in ALL subfolders when searching for tiffs
      'data_path': [current_path], # a list of folders with tiffs 
                                             # (or folder of folders with tiffs if look_one_level_down is True, or subfolders is not empty)                                            
      'subfolders': [], # choose subfolders of 'data_path' to look in (optional)
      'fast_disk': current_path, # string which specifies where the binary file will be stored (should be an SSD)
      'save_folder': current_path
    }

# run one experiment
opsEnd = run_s2p(ops=ops, db=db)