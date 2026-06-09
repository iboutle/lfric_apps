.. -----------------------------------------------------------------------------
    (c) Crown copyright 2026 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _how_to_guide:

How to Guide
===============

This guide provides examples of how to do certain common things in the development of lfric_apps.

.. _create_prog:

Adding a new prognostic variable
----------------------------------------------------------------

Prognostic variables are typically those variables which are required to be passed from one model timestep to the next. However, in lfric_apps, we have a slightly extended definition to include those variables which do not necessarily need to persist from one timestep to the next, but do need to persist for a significant portion of the timestep (e.g. they are passed from slow to fast physics). These fields all live in field collections. Dynamical core variables are created in `create_gungho_prognostics <https://github.com/MetOffice/lfric_apps/blob/main/science/gungho/source/driver/create_gungho_prognostics_mod.F90>`_, and physics variables are created in `create_physics_prognostics <https://github.com/MetOffice/lfric_apps/blob/main/science/gungho/source/driver/create_physics_prognostics_mod.F90>`_.

Adding a new field requires the following code ::

 call processor%apply(make_spec('name', main%field_collection, function_space, &
                                 mult='surface_tiles', twod, is_int,           &
                                 adv_coll=if_adv(advection_flag, adv%last_adv),&
                                 ckp=checkpoint_flag,                          &
				 is_empty=(.not. usage_flag) ) )

where:

* ``name`` is the name you want the new variable to have - it is a text string, e.g. ``tile_temperature``.
* ``field_collection`` is the name of the field collection you want to add the field to, e.g. ``soil_fields``.
* ``function_space`` is the name of the space defining the horizontal and vertical grid the variable is on, e.g. ``WTheta`` or ``W3``.

  * This is optional, and not required if the field is going to be checkpointed, as the information can be inherited from the xml metadata.

* ``mult`` is the name of the multidata dimension you want the field to have.

  * It is an optional input, and not required if the field is going to be checkpointed, as the information can be inherited from the xml metadata. It will also default to just a single dimension field if not provided.

* ``twod`` is a logical flag stating whether the field is 2D (true) or 3D (false).

  * It is an optional input, and not required if the field is going to be checkpointed, as the information can be inherited from the xml metadata. It will also default to false if not provided.

* ``is_int`` is a logical flag stating whether the field is an integer (true) or real (false) field.

  * It is an optional input and will default to false if not provided.

* ``advection_flag`` is a logical flag stating whether a field should be advected or not.

  * It is an optional input, and will default to false if not provided. If it is advected, we also need to specify the collection in which it is advected, e.g. ``adv%last_adv`` to advect on the final outer iteration using the advective form of the transport equation.

* ``checkpoint_flag`` is a logical flag determining whether this field needs adding to the checkpoint-restart dumps, i.e. it needs passing from one timestep to the next.

  * It is an optional input, and will default to false if not provided.

* ``is_empty`` is a logical flag stating whether to create an empty field (true) or allocate memory space for the field (false). If the field is only require under certain logical options, then it should be set to empty when these options are not switched on to save memory space.

  * It is an optional input and will default to false if not provided.

If the field is to be checkpoint-restarted, i.e. you have set ``checkpoint_restart_flag=.true.``, the `lfric_dictionary <https://github.com/MetOffice/lfric_apps/blob/main/applications/lfric_atm/metadata/lfric_dictionary.xml>`_ file must be updated to include the new field which requires checkpointing. The field ``id`` must match the ``name`` argument given above.

Adding a new ancillary field
----------------------------------------------------------------

If you want your new ancillary field to be available throughout the code, then you will need to follow :ref:`these instructions <create_prog>` first to create a prognostic variable to hold the ancillary. However, if you only need the ancillary field to be available for a short time after startup (i.e. it is immediately processed into another field), then this is not required and the field will be made automatically by the steps below.

Ancillary fields are set up in `init_ancils_mod <https://github.com/MetOffice/lfric_apps/blob/main/science/gungho/source/driver/init_ancils_mod.f90>`_ using ::

  call setup_ancil_field( name, depository, ancil_fields, mesh_id, &
                          twod_mesh_id, twod, ndata, time_axis )

where:

* ``name`` you want the new field to have. If this field does not already exist, it will be created. However, if you already created a new prognostic with the same name to hold the ancillary, it will be written to that.
* ``depository`` is the collection which contains all of the fields.
* ``ancil_fields`` is the collection which contains all of the ancillaries to be read.
* ``mesh_id`` and ``twod_mesh_id`` specify the 3D and 2D meshed being used.
* ``twod`` is a logical flag stating whether the field is 2D (true) or 3D (false).

  * It is an optional input which defaults to false.

* ``ndata`` is the size of the non-spatial dimensionality of the field (e.g. number of tiles).

  * It is an optional input which defaults to 1.

* ``time_axis`` is the time axis object for time-varying ancillaries.

  * It is an optional input, the default for which assumes a non-time-varying ancillary.

If you are trying to read a new ancillary field from a pre-existing ancillary file, then the new field must also be added to the appropriate `file_def_ancil <https://github.com/MetOffice/lfric_apps/tree/main/rose-stem/app/lfric_atm/file>`_ file. The field ``id`` is the ``name`` specified above, and the field ``name`` is the name of the field in the netcdf file. N.B. these files live in the workflow to enable flexibility of the input files. The link above points to rose-stem, but for standalone workflows, the lfric app will need modifying in a similar manner.

If you are trying to read an entirely new ancillary file, then several further pieces of infrastructure will be required for this:

1. A new ``<file>`` definition needs adding in the `file_def_ancil <https://github.com/MetOffice/lfric_apps/tree/main/rose-stem/app/lfric_atm/file>`_ files to specify the file to be read.
2. This new file needs adding to `gungho_setup_io <https://github.com/MetOffice/lfric_apps/blob/main/science/gungho/source/driver/gungho_setup_io_mod.F90>`_, where the ``xios_id`` should match the file ``id`` specified in step 1.
3. The new file needs adding to the `rose-meta <https://github.com/MetOffice/lfric_apps/blob/main/science/gungho/rose-meta/lfric-gungho/HEAD/rose-meta.conf>`_ file, the value here being what is used in step 2 to determine where the file is read from.
4. An `upgrade macro <https://github.com/MetOffice/lfric_apps/blob/main/science/gungho/rose-meta/lfric-gungho/versions.py>`_ will be required to add the new file path from step 3 to any configuration namelists (the rose-app.conf).

Adding a new diagnostic
----------------------------------------------------------------

Outputting a prognostic field as a diagnostic should only require one line of code and one line of metadata information. In the code we add ::

  call progfield%write_field('section__name')

where ``progfield`` is the field to be written, ``section`` is the science section containing the diagnostic, and ``name`` is the diagnostic name.

We then must edit the `field_def_diags <https://github.com/MetOffice/lfric_apps/blob/main/applications/lfric_atm/metadata/field_def_diags.xml>`_ to include a matching line for the output diagnostic. The field ``id`` should be the ``section__name`` specified in the code.

However, sometimes we want diagnostics which do not form part of the prognostic evolution of the model. It is important that memory is only allocated for these, and calculations are only done, on the timesteps the diagnostic is requested.
