"""
Example of multiprocessing via a queue based on functions

setproctitle is used for the proces title
"""

import json
import multiprocessing
import os
import time

import setproctitle
from ballistic_integrator.functions.functions import json_encoder


def multiprocessing_job_worker(job_queue, worker_ID, configuration, target_function):
    """
    Function that handles running the job
    """

    # Set name worker
    setproctitle.setproctitle(
        "Disk thickness calculation worker process {}".format(worker_ID)
    )

    # Get items from the job_queue
    for job_dict in iter(job_queue.get, "STOP"):
        #########
        # Stopping or working
        if job_dict == "STOP":
            return None

        ##########
        # Call target function
        print(
            "Worker {} handling system {} {}".format(
                worker_ID,
                "mass_accretor: {}".format(job_dict['mass_accretor']),
                "massratio_accretor_donor: {}".format(job_dict['massratio_accretor_donor']),
            )
        )

        #
        target_function(job_dict)

    return None


def multiprocessing_queue_filler(job_queue, num_cores, configuration, job_iterable):
    """
    Function to handle filling the queue for the multiprocessing
    """

    # Continuously fill the queue
    for job_i, job_dict in enumerate(job_iterable):

        # Put job in queue
        job_queue.put(job_dict)

    # Signal stop to workers
    for _ in range(num_cores):
        job_queue.put("STOP")


def multiprocessing_routine(configuration, job_iterable, target_function):
    """
    Main process to handle the multiprocessing
    """

    print("Starting multiprocessing of grid with the following configuration:")
    print(json.dumps(configuration, indent=4, default=json_encoder))

    ###################
    # Run the calculation through multiprocessing
    start = time.time()

    # Set process name
    setproctitle.setproctitle("Disk thickness calculation parent process")

    # Set up the manager object that can share info between processes NOTE: you can add a result queue here
    manager = multiprocessing.Manager()
    job_queue = manager.Queue(configuration["max_job_queue_size"])

    # Create process instances
    processes = []
    for worker_ID in range(configuration["num_cores"]):
        processes.append(
            multiprocessing.Process(
                target=multiprocessing_job_worker,
                args=(job_queue, worker_ID, configuration, target_function),
            )
        )

    # Activate the processes
    for p in processes:
        p.start()

    # Set up the system_queue
    multiprocessing_queue_filler(
        job_queue=job_queue,
        num_cores=configuration["num_cores"],
        configuration=configuration,
        job_iterable=job_iterable,
    )

    # Join the processes
    for p in processes:
        p.join()
    print(
        "Finished multiprocessing all the tasks. This took {}".format(
            time.time() - start
        )
    )


if __name__ == "__main__":

    configuration = {
        "output_dir": os.path.join(os.getcwd(), "results"),
        "num_cores": 4,
        "max_job_queue_size": 10,
    }

    os.makedirs(configuration["output_dir"], exist_ok=True)

    multiprocessing_routine(configuration)
