.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _jules:

The standalone land surface model - jules
=========================================

The `jules` application runs the JULES land surface model inside LFRic,
without the rest of the atmosphere model. It is meant to reproduce a
standalone JULES point (site) run. Observed near-surface meteorology drives
JULES, and the results can be checked directly against standalone JULES. It
also gives a cheap way to test the LFRic JULES interface.

The application uses the Gungho driver with ``method='jules'`` in the
``formulation`` namelist. Each timestep runs only the JULES surface
calculation (``jules_timestep_alg``). The dynamical core, transport, boundary
layer and radiation schemes are not run. The model is a single column (a
1x1 planar mesh) with one atmospheric layer. That layer represents the
forcing level.

The site forcing data is read through XIOS, and the code that applies it is
in ``interfaces/jules_interface``. Other applications could use it too.

Site forcing
------------

Switching it on
~~~~~~~~~~~~~~~

Two namelists switch the forcing on:

.. code-block:: fortran

   &files
     jules_forcing_path = '/path/to/forcing_file'   ! no .nc suffix
   /
   &initialization
     jules_forcing       = .true.
     jules_forcing_z1_tq = 10.0   ! height of the T and q data (m)
     jules_forcing_z1_uv = 10.0   ! height of the wind data (m)
   /

The forcing file
~~~~~~~~~~~~~~~~

The forcing file is a NetCDF file with one point. By default it uses the
dimensions ``(time, land)`` and must contain the following variables:

=========  ==================================  ==============
Variable   Description                         Units
=========  ==================================  ==============
``Wind``   Wind speed                          m s\ :sup:`-1`
``Tair``   Air temperature                     K
``Qair``   Specific humidity                   kg kg\ :sup:`-1`
``PSurf``  Surface pressure                    Pa
``LWdown`` Downward longwave flux              W m\ :sup:`-2`
``SWdown`` Downward shortwave flux             W m\ :sup:`-2`
``Rainf``  Total precipitation rate            kg m\ :sup:`-2` s\ :sup:`-1`
=========  ==================================  ==============

It must also have a CF ``time`` coordinate variable with units of the form
``seconds since <date>``. The length of the ``time`` dimension sets the size
of the XIOS time axis, and is read when the model starts
(``get_jules_forcing_dim`` in ``time_dimensions_mod``). The ``time``
dimension must be fixed, not unlimited: XIOS treats an unlimited dimension
as its own record dimension rather than as the time axis of the grid, and
fails with "incorrect dimension". A file with an unlimited time dimension
can be converted with ``nccopy -u``. The file does not need any UGRID mesh
information.

The variable names and dimensions are set by the XIOS file definition, so
other formats can be read without converting them. For example,
``file_def_jules_forcing_plumber.xml`` in the rose-stem app reads PLUMBER2
(FLUXNET) met files without renaming or reordering anything: their
variables are on ``(time, y, x)``, and surface pressure and precipitation
are called ``Psurf`` and ``Precip``. Only the unlimited time dimension needs
changing, with ``nccopy -u``.

The XIOS file definition is in ``file_def_jules_forcing.xml``, included by
``iodef.xml``. In the rose-stem app this file is installed from
``file_def_jules_forcing_cardington.xml``, or from
``file_def_jules_forcing_plumber.xml`` by the AR-SLu optional
configuration; the canned example includes
``file_def_jules_forcing_cardington.xml`` directly. It reads the variables
onto an axis-only grid (``jules_forcing_grid``), which avoids needing a
mesh:

.. code-block:: xml

   <file id="jules_forcing" name="jules_forcing" mode="read"
         output_freq="30mi" convention="CF" cyclic="true" enabled=".FALSE.">
     <field id="jules_forcing_time" name="time" unit="seconds"
            operation="once" axis_ref="jules_forcing_time_axis" />
     <field id="forcing_wind_data" name="Wind" operation="instant"
            grid_ref="jules_forcing_grid" />
     ...
   </file>

The file is registered with XIOS in ``gungho_setup_io_mod`` when
``jules_forcing`` is true. Each XIOS field id is the model field name with
``_data`` added, for example ``forcing_wind`` and ``forcing_wind_data``. The
model fields are held in the ``jules_forcing_fields`` collection, which
``create_jules_forcing_fields`` in ``init_jules_forcing_mod`` creates.

.. important::

   The ``output_freq`` of the ``jules_forcing`` file must equal the data
   interval of the forcing file. It must also equal the model timestep.
   The data is not interpolated in time.

Site ancillaries
~~~~~~~~~~~~~~~~

Time-varying leaf area index, as prescribed in standalone JULES, is read
through the ancillary system with ``ancil_option='site'`` in the
``initialization`` namelist. This reads only the leaf area index from the
plant functional type ancillary, ``plant_func_ancil_path`` (relative to
``ancil_directory``, or an absolute path), and updates it as the run
progresses. Unlike the ``fixed`` and ``updating`` options, the data is not
interpolated in time: each value is held until the next data time, as for
the ``'nf'`` option of standalone JULES (``create_fd_ancils`` in
``init_ancils_mod``). The canopy height is not read, and stays at its
``ideal_surface`` value, as in standalone JULES without TRIFFID.

As for the forcing, the site ancillaries are updated before each timestep
(in ``gungho_driver_mod``), so that the data for a time applies from the
start of that timestep. As for prescribed vegetation in standalone JULES,
the leaf area index is then set to at least 0.01 (in
``jules_timestep_alg``). This matters with ``all_tiles`` on, when the
canopy radiation is also calculated for plant types with no leaf area.

The ancillary file does not need any mesh information. In ``iodef.xml``
the ``plant_func_ancil`` file reads the variable ``lai_pft`` on
``(pft, time, x, y)`` onto the axis-only grid ``site_pft_grid``. The
``time`` dimension sets the size of the XIOS axis ``site_ancil_time_axis``.
As for the forcing, the file's ``output_freq`` must equal the timestep.

A standalone JULES prescribed leaf area index file is close to this, but
has its dimensions in the order ``(time, pft, x, y)`` and no time values.
It can be converted with the NCO and netCDF tools, setting the times from
the ``data_start`` and ``data_period`` of the standalone run:

.. code-block:: bash

   ncap2 -O -s 'time=array(0.0,1800.0,$time);
                time@units="seconds since 2010-01-01 00:00:00";
                time@calendar="gregorian"' lai.nc lai_time.nc
   ncpdq -O -a pft,time lai_time.nc lai_pft_time.nc
   nccopy -u lai_pft_time.nc lai_lfric.nc

``ncpdq`` puts the plant type outermost, as LFRic expects for multi-valued
time-varying data, and ``nccopy -u`` makes all the dimensions fixed.

If the file only has records at the times when the leaf area index
changes, it needs a final record after the last timestep of the run.
Otherwise, once the last data time has passed, the time axis wraps round to
the first record and cycles through the data.

Timing
~~~~~~

The forcing time axis is set up without interpolation and with a window size
of one. Each data time is the start of the period that the data applies to,
and the value is held until the next data time. This is the ``'nf'``
interpolation option in standalone JULES. The forcing fields are updated in
``gungho_driver_mod`` after the clock ticks and before ``gungho_step``. As a
result, the first record in the file is used for the first timestep.

How the forcing drives JULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each timestep, before the surface calculations, ``apply_jules_forcing_alg``
calls ``jules_forcing_kernel`` to set the model state from the forcing:

* The surface Exner pressure is calculated from ``PSurf``. The Exner pressure
  at the forcing level comes from hydrostatic balance at height
  ``jules_forcing_z1_tq``.
* Potential temperature in the first layer is ``Tair`` divided by the Exner
  pressure there. The vapour mixing ratio is set from the specific humidity
  ``Qair``, as ``Qair / (1 - Qair)``, because the physics uses mixing ratios.
* The wind ``u`` is set to ``Wind`` and ``v`` is set to zero. Density comes
  from the equation of state.
* The downward SW and LW surface fluxes come from ``SWdown`` and ``LWdown``.
  The photosynthetically active radiation (PAR) is ``visible_fraction``
  (0.5) of the SW. The direct part of PAR is ``1 - diffuse_fraction`` (0.6),
  which matches ``diff_frac = 0.4`` in the JULES kernels.
* ``Rainf`` falls as rain when ``Tair`` is at least 274 K and as snow
  otherwise, as standalone JULES does.

With forcing on, the JULES explicit and implicit kernels use the namelist
forcing heights for ``z1_tq`` and ``z1_uv``, and set ``z1_tq_top`` and
``z1_uv_top`` to zero. ``jules_physics_init`` sets ``i_modiscopt`` to 0,
because the forcing is at a fixed height rather than a model level.

Radiation
~~~~~~~~~

No radiation scheme is run. Instead, ``jules_forcing_init`` (called from
``gungho_model_mod``) sets the SW bands normally read by ``socrates_init``.
It sets two bands, covering the visible (0.2 to 0.69 micron) and
near-infrared (0.69 to 5 micron) parts of the spectrum. The first band has
a ``sw_weight_blue`` of 1 and the second 0. The land albedos do not depend
on the band limits, but the spectral sea albedo (``l_spec_sea_alb``) needs
them: it cannot use more than one band with the broadband limits of -1. With this setup, the existing
``sw_rad_tile`` call in ``rad_tile_alg`` gives the visible direct and
diffuse tile albedos in the first band, and the near-infrared ones in the
second. After ``rad_tile_alg``, ``jules_forcing_sw_up_kernel`` sets the
upward SW flux on each tile to ``SWdown`` multiplied by a weighted sum of the
four albedos, as standalone JULES does. The weights are set by
``wght_alb`` in the ``jules_radiation`` namelist, in the standalone order:
visible direct, visible diffuse, near-infrared direct, near-infrared
diffuse. The standalone default of ``0.0,0.5,0.0,0.5`` uses only the diffuse
albedos, weighted equally.

The downward LW comes straight from the forcing. The tile emissivities come
from the existing LW tile calculation.

The solar zenith angle comes from ``illuminate_alg``, which is called in
``jules_timestep_alg`` each timestep. The radiation, orbit and star namelists
are not needed: ``illuminate`` uses default orbital parameters and a
radiation step of one when the radiation scheme is not Socrates.

Initial state and surface properties
------------------------------------

There are no ancillary files. The surface is set from the ``ideal_surface``
namelist, and the soil from the analytic branch of ``init_soil_fields_alg``.
To support site runs, these ``ideal_surface`` items were added. They were
previously hard-coded:

* soil hydraulic and thermal properties: ``clapp_horn_b``, ``soil_albedo``,
  ``soil_cond_sat``, ``soil_moist_crit``, ``soil_moist_sat``,
  ``soil_moist_wilt``, ``soil_suction_sat``, ``soil_thermal_cap`` and
  ``soil_thermal_cond``;
* the TOPMODEL topographic index mean and standard deviation,
  ``mean_topog_index`` and ``stdev_topog_index``;
* the initial ``water_table`` and ``wetness_under_soil``.

The ``vn32`` upgrade macros in ``jules-lsm`` and ``lfric-gungho`` set the
new items to the old hard-coded values. This leaves existing configurations
unchanged. The same applies to these JULES settings, which were hard-coded
before and are now in the namelists of all applications using JULES:

* ``lai_alb_lim_sn`` in ``jules_snow``, the minimum leaf area index of
  snow-covered canopies, with one value for each plant type;
* ``soilhc_method`` in ``jules_soil``, the soil thermal conductivity model
  (``'cox'``, ``'johansen'`` or ``'chadburn'``);
* ``wght_alb`` in ``jules_radiation``, only used with site forcing.

Plant functional types
~~~~~~~~~~~~~~~~~~~~~~

As well as the 5 standard plant types (``brd_leaf``, ``ndl_leaf``,
``c3_grass``, ``c4_grass`` and ``shrub``), the ``jules_surface_types``
namelist accepts the 13 plant types of standalone JULES:
``brd_leaf_dec``, ``brd_leaf_eg_trop``, ``brd_leaf_eg_temp``,
``ndl_leaf_dec``, ``ndl_leaf_eg``, ``c3_grass``, ``c3_crop``,
``c3_pasture``, ``c4_grass``, ``c4_crop``, ``c4_pasture``, ``shrub_dec``
and ``shrub_eg``. Each plant type has a ``jules_pftparm`` namelist,
identified by ``pft_name_io``. Whether snow can be held under the canopy is
set by ``cansnowpft_io`` in each of these, rather than by the ``cansnowpft``
array of standalone JULES's ``jules_snow`` namelist, which LFRic could only
read with 5 plant types. The metadata for the extra types is in
``jules-lsm``; it could move to the shared JULES metadata in due course.
Surface types that are not used are left unset.

The atmospheric profile namelists only use ``profile_size=1`` placeholder
values, because the forcing overwrites the first layer before JULES runs.
``map_physics_fields_alg`` is skipped when forcing is on.

Configuration
-------------

Required namelists
~~~~~~~~~~~~~~~~~~

``applications/jules/source/jules_mod.f90`` defines
``jules_required_namelists``, which ``jules.f90`` uses in place of the Gungho
list. Only these namelists are required: ``finite_element``,
``formulation``, ``base_mesh``, ``initial_wind``, ``planet``,
``timestepping``, ``extrusion``, ``orography`` and ``energy_correction``.

The ``energy_correction`` namelist is required because a missing
``encorr_usage`` would switch energy correction on. That then fails on a
planar mesh.

The ``lfric-jules`` metadata imports ``lfric-gungho`` and makes the unused
namelists optional (``compulsory=false``). These are ``boundaries``,
``checks``, ``convection``, ``cosp``, ``departure_points``, ``esm_couple``,
``helmholtz_solver``, ``iau``, ``initial_density``, ``mixed_solver``,
``mixing``, ``solver``, ``specified_surface`` and ``transport``. The
transport checks in ``check_configuration_mod`` only run when the
``transport`` namelist has been loaded. The model is run with
``boundary_layer='none'``.

The canned example
~~~~~~~~~~~~~~~~~~

``applications/jules/example`` has a ``configuration.nml`` and
``iodef.xml`` for the Cardington site (52.103 N, 359.576 E). The setup is:

* start time 2005-01-01 00:00;
* 143 timesteps of 1800 s;
* 100% C3 grass;
* 10 m forcing heights.

The initial state, soil properties and topographic index match the
standalone JULES suite u-dz368 (RAL3). Before running, set
``jules_forcing_path`` to the location of the Cardington forcing file.

Rose stem
~~~~~~~~~

The rose-stem app is ``rose-stem/app/jules``. The Cardington test runs as
the task ``cardington-1x1-50000x50000`` in ``tasks_jules.cylc``. It is
included in the group
``jules_cardington-1x1-50000x50000_azspice_gnu_fast-debug-64bit``. The main
app is kept consistent with ``lfric_atm``. For Cardington, three optional
configurations apply, in this order:

``rose-app-1x1-50000x50000.conf``
   Mesh and ``base_mesh`` settings.

``rose-app-ral4.conf``
   The JULES settings from the ``lfric_atm`` RAL3 optional configuration
   (``rose-stem/app/lfric_atm/opt/rose-app-ral3.conf``). It covers the lake,
   plant functional type and snow parameters, plus the ``jules_surface``
   switches. The ``lfric_atm`` settings for other schemes are left out,
   because those schemes are not run in this app.

``rose-app-cardington.conf``
   The settings where u-dz368 differs from ``ral4``.
   It is applied after ``ral4`` and overrides it. The settings are:

   * plant functional type albedo parameters (``alnir``, ``alpar``,
     ``omega``, ``omnir``);
   * shrub roughness (``z0hm_pft_io=0.05``, ``z0v_io=0.2``), which matches
     the RAL3 configuration of u-dz368;
   * snow canopy parameters (``can_clump`` and ``n_lai_exposed``);
   * ``formdrag='none'``;
   * ``l_mo_buoyancy_calc``, ``l_var_rainfrac`` and ``l_sea_alb_var_chl``
     switched off;
   * a CO\ :sub:`2` mixing ratio of 5.941e-4.

A second test replicates the PLUMBER2 site AR-SLu (San Luis, Argentina) from
the standalone JULES suite u-dy575. It runs for the whole of 2010 as the
task ``ar_slu-1x1-50000x50000``, with the mesh optional configuration and
``rose-app-ar_slu.conf``. That optional configuration sets:

* the 13 plant types with the u-dy575 parameters, and 5 non-vegetated types
  including the urban canyon and roof (``l_urban2t``), giving 18 land tiles;
* the site tile fractions, soil properties and topographic index;
* the initial state, taken from the u-dy575 dump at the start of its main
  run, after spin-up (soil moisture and temperature, water table, wetness
  below the soil column and tile temperatures);
* time-varying leaf area index from a site ancillary
  (``ancil_option='site'``), the MODIS leaf area index file used by
  u-dy575 converted as described above;
* the u-dy575 canopy heights, which are constant;
* forcing heights of 7.985 m, the u-dy575 albedo weights
  (``0.3,0.2,0.34,0.16``) and the Chadburn soil thermal conductivity;
* ``f_lon_deg=0.0``. The forcing is in local time, and u-dy575 uses
  ``l_local_solar_time``, so the solar hour angle must not include the
  longitude.

The site has no urban fraction, so the MORUSES urban options are left off;
they would otherwise need urban ancillary fields. As in u-dy575, the
rainfall is treated as point data (``l_point_data``), which changes the
canopy interception and surface runoff.

The PLUMBER2 met file is read using ``file_def_jules_forcing_plumber.xml``,
after making its time dimension fixed with ``nccopy -u``.

The ``[file:configuration.nml]`` source list holds every namelist, including
the ones that are ignored. Optional ones are marked ``(optional)``. This means
a namelist that is trigger-ignored is still output if a setting brings it
back into use.

Comparison with standalone JULES
--------------------------------

The Cardington test has been compared with the equivalent standalone JULES
run (u-dz368) for January 2005. The two agree closely:

* soil moisture, soil temperature and runoff (including TOPMODEL baseflow)
  are effectively identical;
* surface energy-balance terms agree to within about 0.3 W m\ :sup:`-2`.

There is one known remaining difference. Soil evaporation is about 2% higher
in LFRic, which offsets the latent and sensible heat fluxes by
+0.3 W m\ :sup:`-2` and -0.2 W m\ :sup:`-2`. This difference is present from
the first timestep, and its cause has not yet been found. The 1.5 m specific
humidity is also slightly lower, by about 0.02 g kg\ :sup:`-1`. This is
expected, because the screen diagnostic method (``iscrntdiag``) cannot be
chosen in LFRic.

Limitations
-----------

* Only single-point runs on one MPI rank are supported.
* The model timestep must equal the interval of the forcing data, and the
  ``output_freq`` of the ``jules_forcing`` file must match it. The data is
  not interpolated in time.
* All precipitation is supplied as ``Rainf`` and split into rain and snow at
  274 K. Separate rain, snow, or convective and large-scale precipitation
  are not supported.
* Point rainfall data (``l_point_data``) is only available with site
  forcing.
* Some standalone JULES options cannot be set in LFRic. These include
  ``iscrntdiag``, tiled soil (``l_tile_soil``), ``l_local_solar_time``
  (see the AR-SLu configuration for the workaround) and the
  ``jules_temp_fixes`` switches, which are hard-coded in
  ``jules_physics_init``.
* There is no spin-up. A spun-up state has to be supplied as the initial
  state in ``ideal_surface``.
* Diffuse fraction and visible fraction are fixed at 0.4 and 0.5.
* Some JULES inputs are deliberately left unset, because the LFRic build of
  JULES does not use them. Examples are ``ti_sig_soilt``, ``a_step``, and the
  momentum terms of ``surf_couple_implicit``.
* Many fields that are only used by the rest of the atmosphere model are not
  initialised. They should not be written as diagnostics.
