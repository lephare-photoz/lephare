import multiprocessing
import os


def _set_omp_num_threads():
    # Check the maximal number of CPUs
    try:
        max_cpus = multiprocessing.cpu_count()
    except NotImplementedError:  # pragma no cover
        return

    # If the environment variable is not defined
    # Take the max between 1 and max_cpus-1
    if "OMP_NUM_THREADS" not in os.environ:
        # Si non défini, utilise max(1, max_cpus - 1)
        omp_num_threads = max(1, max_cpus - 1)
        os.environ["OMP_NUM_THREADS"] = str(omp_num_threads)
    else:
        # If the environment variable is defined
        # Take the max between 1 and max_cpus-1 if out the limits
        try:
            user_threads = int(os.environ["OMP_NUM_THREADS"])
            if user_threads > max_cpus:
                omp_num_threads = max(1, max_cpus - 1)
                os.environ["OMP_NUM_THREADS"] = str(omp_num_threads)
            elif user_threads < 1:
                omp_num_threads = 1
                os.environ["OMP_NUM_THREADS"] = str(omp_num_threads)
            else:
                omp_num_threads = user_threads
        except ValueError:
            # If not a valid integer
            omp_num_threads = max(1, max_cpus - 1)
            os.environ["OMP_NUM_THREADS"] = str(omp_num_threads)

    print(f"OMP_NUM_THREADS set to {os.environ['OMP_NUM_THREADS']} (max available CPUs: {max_cpus})")
