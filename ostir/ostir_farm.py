# This module contains the rewrite of the ostir factory. It enables running multiple sequences at one.
# This has the result of cleaning up the CLI for the progress bar.


# ------------ Imports ------------
import importlib.util
import concurrent.futures
from typing import Callable
from dataclasses import dataclass

from .ostir_calculations import find_start_codons
from .ostir_worker import ostir_worker
from .data_classes import ostir_task


# Optional dependency
if importlib.util.find_spec("progress") is not None:
    from progress.bar import Bar
else:
    Bar = None

# ------------ Supporting Code ------------


def _ostir_progress_callback(tasks: int) -> Callable:
    # Update the progress bar
    if Bar is not None:
        bar = Bar('Running OSTIR RBS Predictions: ', max=tasks)
    else:
        return None
    max = tasks
    itterations = 0

    def _callback():
        nonlocal itterations
        itterations += 1
        nonlocal max
        if itterations == max:
            bar.next()
            bar.finish()
            return None

        bar.next()
        return None


    return _callback


# -------------- Main Code ------------
def ostir_farm(task_objects: list, threads: int=1, callback: Callable=None, verbosity: int=0) -> dict:
    """ Run OSTIR on multiple sequences in parallel. """

    # Find all start codons in every supplied sequence.
    # Note: returns list of generators
    start_codons = [find_start_codons(task.sequence, start_range=[task.start, task.end]) for task in task_objects]
    id_values = [i for i in range(len(task_objects))]

    # Tie the generators to the task objects
    sets = [(task, start_codon, id_value) for task, start_codon, id_value in zip(task_objects, start_codons, id_values)]

    # Unpack the generators
    arguments = [(task, start_codon, id_value, verbosity) for
                  task, start_codons, id_value in sets for start_codon in start_codons]



    # Tell the callback function how many tasks we have. Expect back another callback function.
    if callback is not None:
        callback = callback(len(arguments))
    elif Bar and verbosity > 0:
        callback = _ostir_progress_callback(len(arguments))
    else:
        callback = lambda: None

    # Create a pool of workers and run OSTIR

    if threads > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as multiprocessor:
            #parallel_output = multiprocessor.map(self._parallel_dG, *parallelizer_arguments)

            futures: list = []
            for task in arguments:
                future_object = multiprocessor.submit(ostir_worker, *task)
                future_object.add_done_callback(lambda x: callback())
                futures.append(future_object)
            parallel_output = [future.result() for future in futures]
    else:
        parallel_output: list = []
        for task in arguments:
            result = ostir_worker(*task)
            parallel_output.append(result)
            callback()

    parallel_output = [x for x in parallel_output if x is not None]

    # Sort out to each ID back the results
    results = {}
    for [id, result] in parallel_output:
        if id not in results:
            results[id] = {'name': result.name,
                                        'results': []}
        results[id]['results'].append(result)

    # Sort the results within each ID
    for item in results.items():
        item[1]['results'].sort(key=lambda x: x.start_position)

    # Return the results
    return parallel_output
    
# Don't forget the end of file newline!
