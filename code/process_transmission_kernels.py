import scipy.stats as stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.interpolate as interp
import moveSTIR as stir
import os

"""
Process transmission matrices
"""

def interpolate_all_trajectories(step, dat):
    """
    Linearly interpolate all of the pig trajectories to the step scale.
    An alternative approach to fitting the CTMM model.

    Parameters
    ----------
    step : int
        Step size on which to linear interpolate pig trajectories
    dat : DataFrame
        All pig data 

    Returns
    -------
    : dict
        Keys are host IDs, values are interpolated movement trajectories
        to the "step" scale.  
    """

    unq_collar = np.unique(dat.individual_ID)

    # Ensure all pigs are aligned when interpolating
    interp_vals = np.arange(dat.unix_time.min(), dat.unix_time.max() + step, step=step)

    all_fitted = {}
    for unq_ind in unq_collar:

        trial_dat = dat[dat.individual_ID == unq_ind]

        # Remove any of the same datetimes
        trial_dat = trial_dat[~trial_dat.datetime.duplicated()].sort_values("datetime").reset_index(drop=True)

        min_time = trial_dat.unix_time.min()
        time = trial_dat.unix_time.values
        xloc = trial_dat.UTMx.values
        yloc = trial_dat.UTMy.values

        fitted = stir.fit_interp_to_movement(time, xloc, yloc, interp_vals=interp_vals)
        all_fitted[unq_ind] = fitted

    return(all_fitted)


def process_transmission_kernels(all_fitted, params, memory_map=True,
                                 path=os.getcwd(), save_append=""):
    """
    Process transmission kernels to extract the marginals (row and column) sums.

    Summing across rows (axis=1) is the deposition marginal 
        - Deposition marginal: How a host's position at time t contributes to FOI 
        in all other hosts (a rate)
    Summing across columns (axis=0) is the foi marginal
        - The force of infection experienced by a host at time t from all 
        other hosts (a rate)

    Parameters
    ----------
    all_fitted : dict
        Keys are host IDs that look up data frames that contain the discretized
        CTMMs for the each host. NOTE: Make sure predicted CTMM values are at
        the same discretized time steps for all hosts being compared.
    params : dict
        beta : float, acquisition rate
        lam : float, deposition rate
        pathogen_decay : float, Path decay rate
        distance_decay : float, Either max distance (cutoff) or mean distance (gaussian)
                         of contact function
        dd_type : str, 'cutoff' or 'gaussian'
    memory_map : bool
        If True use memory mapping to process transmission kernels.  This is
        necessary when trajectories are longer than about 20,000 steps as 
        the transmission matrix cannot be held in RAM.  This function will 
        save the transmission kernel to disk and use memory mapping to
        manipulate it. 

        WARNING: These matrices can be LARGE (potentially 
        hundreds of GB or larger). Make sure you have
        sufficient space on your disk or on an external hard drive. The 
        transmission matrices are removed after processing the marginals.
        If False, manipulates the transmission kernel in memory.  
    path : str
        The path where the temporary memory mapped transmission kernel should
        be saved, as well as the resulting marginal calculations.  Default
        is the current working directory.
    save_append : str
        A string that will be appended to the saved results.

    Returns
    -------
    None
        Pickled results are saved to disk as 
        `marginal_fois_{0}.pkl`.format(save_append).  The pickled results
        are a tuple that contain
            (deposition marginal, foi direct marginal, foi indirect marginal,
             host keys, parameter dictionary)
    """

    host_keys = list(all_fitted.keys())
    n = len(host_keys)

    beta = params['beta']
    lam = params['lam'] 
    pathogen_decay = params['pathogen_decay'] 
    distance_decay = params['distance_decay'] 
    dd_type = params['dd_type'] 

    # Loop through pairwise sums
    count = 1
    fois_direct = {}
    fois_indirect = {}
    deposit = {}

    for i, h1_nm in enumerate(host_keys):
        for j, h2_nm in enumerate(host_keys):
	    
            if i != j:
                print("{0} of {1}".format(count, n*(n - 1)))
                print("Working on {0} lag {1}".format(h1_nm, h2_nm))
                host1 = all_fitted[h1_nm].copy()
                host2 = all_fitted[h2_nm].copy()

                # Align time stamps. We are only comparing hosts where they overlap in time
                mintime = np.max([np.min(host1.time), np.min(host2.time)])
                maxtime = np.min([np.max(host1.time), np.max(host2.time)])
                host1 = host1[(host1.time >= mintime) & (host1.time <= maxtime)].reset_index(drop=False)
                host2 = host2[(host2.time >= mintime) & (host2.time <= maxtime)].reset_index(drop=False)
                nt = len(host1)

                try:
                    if memory_map:
                        # Use memory mapping for matrices that can't be stored
                        # in memory
                        datapath = os.path.join(path, "K_{0}lag{1}.dat".format(h1_nm, h2_nm))
                        deltat, K_1lag2, _ = stir.transmission_kernel(host1, host2, pathogen_decay, distance_decay, 
                                                                      dd_type=dd_type, max_size=1,
                                                                      file_path=datapath)

                    else:
                        deltat, K_1lag2, _ = stir.transmission_kernel(host1, host2, pathogen_decay, distance_decay, 
                                                                      dd_type=dd_type, max_size=20000)

                    tres = summarize_K(K_1lag2, deltat, nt, beta, lam) 

                    # Remove Kmat from disk
                    if memory_map:
                        try:
                            os.remove(datapath)
                        except FileNoteFoundError:
                            pass

                    time_range = host1.time.max() - host1.time.min()
                    deposit[(i, j)] = (tres[0], time_range, deltat, host1.time.values)
                    fois_direct[(i, j)] = (tres[1], time_range, deltat, host1.time.values)
                    fois_indirect[(i, j)] = (tres[2], time_range, deltat, host1.time.values)

                except IndexError:
                    print("Error for {0} and {1}. No overlap".format(h1_nm, h2_nm))

                count += 1

    # Save results
    pd.to_pickle((deposit, fois_direct, fois_indirect, host_keys, params), 
                 os.path.join(path, "marginal_fois_{0}.pkl".format(save_append)))
    return(None)


def summarize_K(K, deltat, n, beta, lam):
    """
    Given a transmission matrix K (in memory or memory mapped), compute
    the marginal summaries

    Parameters
    ----------
    K : array or str
    	If array, this is the transmission kernel
    	If str, this is the file path to the transmission kernel on disk
    deltat : float
    	Time step
    n : float
    	Size of K matrix
    beta : float
    	Transmission parameter
    lam : float
    	Shedding parameter

    Return
    ------
    : tuple
    	(deposition marginal, foi direct marginal, foi indirect marginal)
    	Direct marginal only extracts the diagonal

    """

    if type(K) == str:
        # Matrix is stored on disk
        Kmat = np.memmap(K, dtype=np.float32, shape=(n, n))
    else:
        # Matrix is an array stored in memory
        Kmat = K

    deposit = beta * lam * (Kmat*deltat).sum(axis=0)
    foi = beta * lam * (Kmat*deltat).sum(axis=1)
    foi_direct = beta * lam * (np.diagonal(Kmat)*deltat)
    foi_indirect = foi - foi_direct

    res = (deposit, foi_direct, foi_indirect)
    return(res)


if __name__ == '__main__':
	
	# Load the pig data
    dat = pd.read_csv("../data/pig_movements.csv")
    dat = (dat.assign(datetime=lambda x: pd.to_datetime(x.date_time))
              .assign(unix_time = lambda x: x.datetime.astype(np.int64) / (60 * 10**9)))

    # Drop host 12...there is no temporal overlap
    dat = dat[dat.individual_ID != 12]
    host_keys = np.unique(dat.individual_ID)

	# Extract/compute interpolated trajectories
    step = 5
    use_ctmm = True

    if use_ctmm:
        all_fitted = {}
        for h in host_keys:
            af = pd.read_csv("../data/ctmm_data/traj_{0}.csv".format(h))
            af = af.assign(time=lambda x: x.time_s / 60)
            af = af[["x", "y", "time"]]
            all_fitted[h] = af
    else:
        all_fitted = interpolate_all_trajectories(step, dat)

    host_keys = list(all_fitted.keys())

	# Process transmission kernels
    params = dict(beta = 1, # Just need a proportion so absolute value doesn't matter here
                  lam = 1, # Just need a proportion so absolute value doesn't matter here
                  pathogen_decay = 1 / (24 * 60. * 5), # Average persistence pathogen of 5 day on the minute scale
                  distance_decay = 10, # meters
                  dd_type = "cutoff")

    # all_fitted_red = {1: all_fitted[1], 27: all_fitted[27]}
    process_transmission_kernels(all_fitted, params, memory_map=True,
                                 path="../results/trans_kernels",
                                 save_append="step{0}_use_ctmm_{1}".format(step, use_ctmm))




