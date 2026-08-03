##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
##############################################################################
# Some of the content of this file has been produced with the assistance of
# Anthropic Claude Opus 5 (Claude Code).
##############################################################################

'''
PSyclone transformation script for the LFRic (Dynamo0p3) API to apply
colouring and redundant computation to the level-1 halo for
the initialisation built-ins generically.

'''

from psyclone_tools import (redundant_computation_setval, colour_loops,
                            view_transformed_schedule)

from psyclone.psyir.transformations import OMPParallelTrans
from psyclone.transformations import LFRicOMPLoopTrans

# The CASIM mechanistic activation keeps its state in module level variables
# inside CASIM, so the invoke containing it must not be threaded over cells.
NO_OMP_KERNEL = "casim_aerosol_act_code"


def trans(psy):
    '''
    Applies PSyclone colouring and redundant computation transformations.

    '''
    redundant_computation_setval(psy)
    colour_loops(psy)

    # Extracted from psyclone tools#
    # openmp_parallelise_loops#
    # To avoid adding OMP around the loop containing the CASIM activation

    otrans = LFRicOMPLoopTrans()
    oregtrans = OMPParallelTrans()

    # Loop over all the Invokes in the PSy object
    for invoke in psy.invokes.invoke_list:
        schedule = invoke.schedule
        if any(kern.name.lower() == NO_OMP_KERNEL
               for kern in schedule.coded_kernels()):
            continue
        # Add OpenMP to loops unless they are over colours or are null
        for loop in schedule.loops():
            if loop.loop_type not in ["colours", "null"]:
                oregtrans.apply(loop)
                otrans.apply(loop, options={"reprod": True})
    # Extracted from psyclone tools#

    view_transformed_schedule(psy)
