import numpy as np
import sys
# option to import from github folder
sys.path.append('C:/Users/JJB/Anaconda3/envs/suite2p/Lib/site-packages/suite2p')
from suite2p import run_s2p
sys.path.append('C:/ASDROOT/STUDY/TwoPhotonDataAnalysis')
#print(sys.path)
from forSuite2p import mydefault_ops
import glob
import os



#current_rootPath = 'C:/ASDROOT/STUDY/twoPhotonData_motionCorrected/113Recording_20230123A_Ding_Site09B_sameFOV0122'
#current_rootPath = 'C:/ASDROOT/STUDY/twoPhotonData_motionCorrected/113Recording_20230122A_Ding_Site09B'
current_rootPath = os.getcwd()
print('current_rootPath is %s\n' %(current_rootPath),end="")

target_dir_prefix = 'Result'


file_list = glob.glob((current_rootPath+'/'+target_dir_prefix+'*'))
print('Find %d dir!\n' %(len(file_list)),end="")
for tempi in range(len(file_list)):
   temp_fullDir = os.path.split(file_list[tempi])
   temp_dir = temp_fullDir[1]
   print('%s\n' %(temp_dir),end="")
   
target_fullDir = os.path.split(file_list[len(file_list)-1])
target_dir = target_fullDir[1]

#target_dir = 'Result20230128T170442'

print('Choose %s to analyse!\n' %(target_dir),end="")
   
   



current_path = current_rootPath + '/' + target_dir

os.chdir('C:/ASDROOT/STUDY/TwoPhotonDataAnalysis/forSuite2p')
current_rootPath = os.getcwd()


# set your options for running
ops = mydefault_ops.default_ops() # populates ops with the default options

ops['threshold_scaling'] = 0.45 # 1.5-->0.75-->0.70-->0.52-->0.45
ops['connected'] = False  # whether or not to keep ROIs fully connected (set to 0 for dendrites)
ops['max_iterations'] = 8 # maximum number of iterations to do cell detection 20-->40-->20-->6-->8
ops['batch_size'] = 500 # 500-->1500-->500
ops['circular_neuropil'] = True # to enable circular extend of neuropil duration neuropil search
ops['active_percentile'] = 99 # to search roi by a percentile based-on distribution, instead of a fixed parameter, 95-->99


#ops['nbinned'] = 2500 # 10000 --> 2500

#ops['anatomical_only'] = 4 # run cellpose to get masks on 1: max_proj / mean_img; 2: mean_img; 3: mean_img enhanced, 4: max_proj
#ops['diameter'] = 0


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
    
    
filename = 'C:\ASDROOT\STUDY\TwoPhotonDataAnalysis\my2pScripts\m_runScriptState.dat'
m_runScriptState = np.memmap(filename, dtype='float64', mode='r+', shape=(1, 1))
m_runScriptState[0,0] = 0

# run one experiment
opsEnd = run_s2p(ops=ops, db=db)


m_runScriptState[0,0] = 1
    
    
print('suite2p over\n',end="")
    
# t = time.time()
# waitTime = 2 # secs

# while True
#    elapsed = time.time() - t
#    if elapsed > waitTime
#       print('suite2p over, additional waitTime = %.2f secs\n' %(waitTime),end="")
#       break




