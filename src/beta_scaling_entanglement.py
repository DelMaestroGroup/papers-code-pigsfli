# entanglement_combiner.py (now also uses jacknife to get S2)

# Takes average <S2> from many random seeds
# and combines them into one file

# Made for files in which SWAP histograms are of size 1

import os
import numpy as np

#  Set value of U
# U_list=np.array([0.500000])
# U_list=np.array([3.300000])
U_list=np.array([10.000000])

beta_list=[0.6,0.7,0.8,0.9,1.0,1.15,1.3,1.5,1.75,2.0,2.5,3.0,3.5,4.0]


# bin size
bs = 10000

# dimension of the hypercube (1,2,or3)
D = 1

# linear size of the hypercube
L = 8

# total number of particles in the system
N = 8

# maximum linear size of the subregion
l_max = 4

# Append sweep results to same list so we can copy paste to plotting script
S2_plot = []
S2_err_plot = []
for U in U_list:
    for beta in beta_list:
        for lA_sector_wanted in [4]:
            
            # total number of sites in the subregion (square geom.)
            mA_sector_wanted = lA_sector_wanted**D
            
            incomplete_seeds = [] 
            seeds_list = list(range(1000))
            seeds_measured = []

            # Set desired parameters of calculation
            L = "%d"%(L)
            N = "%d"%(N)
            l_max = "%d"%l_max
            beta = "%.6f"%beta
            bin_size = "%d"%bs
            D = "%d"%D
            U = "%.6f"%(U)
            t = "1.000000"

            # Get path where raw data for the simulation is stored
            path = "/Users/ecasiano/Desktop/PaperData/PaperData/"
            path += D+"D_"+L+"_"+N+"_"+l_max+"_"+U+"_"+\
            t+"_"+beta+"_"+bin_size+"/"

            # Stores all files in the directory
            filenames_all = os.listdir(path)

            # Convert strings back to numbers
            L = int(L)
            N = int(N)
            l_max = int(l_max)
            beta = float(beta)
            bin_size = int(bin_size)
            D = int(D)
            U = float(U)
            t = float(t)

            # Saves the files relevant to the Renyi Entanglement Entropy calculation
            files_SWAP = []

            # Iterate over all filenames
            for filename in filenames_all:

                # Extract parameter information from file name
                parameters = filename.split("_")

                if parameters[0]=='1D' or parameters[0]=='2D' or parameters[0]=='3D':

                    D_want = int((parameters[0])[0]) # hypercube dimension
                    L_want = int(parameters[1]) # hypercube linear size
                    N_want = int(parameters[2]) # total particles
                    l_want = int(parameters[3]) # subsystem linear size (actually l_max)
                    U_want = float(parameters[4]) # interaction potential
                    t_want = float(parameters[5]) # tunneling parameter
                    beta_want = float(parameters[6]) # imaginary time length (K_B*T)**(-1)
                    bin_size_want = int(parameters[7])
                    filetype = (parameters[8]) # identifies the data stored in file
                    seed = int(parameters[9].split(".")[0]) # random seed used

                    if filetype=='SWAP':

                        # Set parameters of simulations from differenet seeds we want to evaluate [D,L,N,l,U,t,beta,bins,type]
                        parameters_to_evaluate = [D_want,
                                                  L_want,
                                                  N_want,
                                                  l_want,
                                                  U_want,
                                                  t_want,
                                                  beta_want,
                                                  bin_size_want,
                                                  'SWAP']

                        if [D,L,N,l_max,U,t,beta,bin_size,filetype] == parameters_to_evaluate:
                            if os.stat(path+filename).st_size > 0:
                                with open(path+filename) as f:
                                   count = sum(1 for _ in f)
                                if count > 1: # only consider files that managed to save at least 100 bins
                                    files_SWAP.append(filename)
                                    seeds_measured.append(seed)
                                else:
                                    incomplete_seeds.append(seed)

            # Sort SWAP files in ascending order of random seed
            files_SWAP.sort()

            # Get total number of seeds and total columns per data file
            reference_file = np.loadtxt(path+files_SWAP[0])
            number_of_seeds = len(files_SWAP)
            columns_per_file = reference_file.shape[1]
            
            print("\nN: ",N)
            print("L: ",L)
            print("l_A: ", lA_sector_wanted)
            print("U: ",U)
            print("beta: ",beta)

            # Get column sum of SWAP files for each seed
            m_max = l_max**D # total number of sites in SWAP region
            SWAP_col_sums = np.zeros((number_of_seeds,m_max+1)).astype(int)
            combined_SWAP_data = np.zeros((number_of_seeds,columns_per_file))
            for i,filename in enumerate(files_SWAP):
                data = np.loadtxt(path+filename)
                data_mean = np.mean(data,axis=0)
                combined_SWAP_data[i] = data_mean 
                
                SWAP_col_sums[i] = np.sum(data,axis=0)
                
            # --- Jacknife --- #
            SWAP_col_seed_sum = np.sum(SWAP_col_sums,axis=0)
            
            # Initialize structure where jacknifed S2 values will be stored
            S2_jacknifed = np.zeros(SWAP_col_sums.shape)
            
            for i in range(SWAP_col_sums.shape[0]):
                
                # Generate jacknifed data
                SWAP_jacknifed_sum = SWAP_col_seed_sum - SWAP_col_sums[i]
                
                SWAP_m = SWAP_jacknifed_sum
                SWAP_0 = SWAP_jacknifed_sum[0]
                S2_jacknifed[i] = -np.log(SWAP_m/SWAP_0)

            # Get mean and std dev,err of S2
            S2_mean = np.mean(S2_jacknifed,axis=0)
            S2_stderr = np.std(S2_jacknifed,axis=0)  * np.sqrt(S2_jacknifed.shape[0])

            print("<S2> = %.4f +/- %.4f"%(S2_mean[mA_sector_wanted],S2_stderr[mA_sector_wanted]))

            # Save the l=2 results
            S2_plot.append(S2_mean[mA_sector_wanted])
            S2_err_plot.append(S2_stderr[mA_sector_wanted])
            
    print("\n\n")
    for result in S2_plot:
        print(f"{result:0.8f}",end=",")
    print("\n\n")
    for error in S2_err_plot:
        print(f"{error:0.8f}",end=",")

    print("\n")
    print("<S2>=",S2_plot)
    print("S2_err=",S2_err_plot)
    print("beta=",beta_list)
    
    print("mA_sector_wanted ",mA_sector_wanted)


print("number of seeds: ",number_of_seeds)
print("incomplete seeds: ",[int(i) for i in incomplete_seeds])
#seeds_measured.sort()
#print([x for x in range(seeds_measured[0],seeds_measured[-1]+1) if x not in seeds_measured])

beta_list = np.array(beta_list)
 #Format the data file
with open("../ProcessedData/"+str(D)+"D_%d_%d_%d_%.6f_%.6f_betas_%d_S2.dat"%(L,N,l_max,U,t,bin_size),"w+") as processed_data:
    np.savetxt(processed_data,np.c_[beta_list,S2_plot,S2_err_plot],delimiter=" ",fmt="%.16f",header="BH Parameters: L=%d,N=%d,D=%d,l=%d,U=%.6f,t=%.6f,bin_size=%d \n beta            <S2>               StdErr."%(L,N,D,l_max,U,t,bin_size))