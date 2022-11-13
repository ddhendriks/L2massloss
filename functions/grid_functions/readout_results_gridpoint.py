def readout_results_gridpoint():
    """
    Function to read out the results of each gridpoint
    """

    filename = "L2massloss/disk_sltns/fL2grid_M1.4_q0.5_case3.pkl"
    with open(filename, 'rb') as f:
        object_file = pickle.load(f)

    #
    header = object_file[0]
    data_list = object_file[1:]
    data_array = np.zeros((data_list[0].shape[0], len(data_list))).T

    for arr_i, arr in enumerate(data_list):
        print(arr.shape)
        # data_array[arr_i][:] = arr[:]

    # df = pd.DataFrame(data_array, columns=header)
