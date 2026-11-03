import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version31_32 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


"""
Copy this template and complete to add your macro

class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>

    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"

    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""


class vn32_tXXX(MacroUpgrade):
    """Upgrade macro for ticket #XXX by Ian Boutle."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_tXXX"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lsm
        """Add the soil properties, topographic index, initial water table and
        wetness under the soil column, with the values previously hard-coded"""
        for name, value in [
            ("clapp_horn_b", "11.2"),
            ("mean_topog_index", "1.0"),
            ("soil_albedo", "0.11"),
            ("soil_cond_sat", "0.0015"),
            ("soil_moist_crit", "0.37"),
            ("soil_moist_sat", "0.456"),
            ("soil_moist_wilt", "0.263"),
            ("soil_suction_sat", "0.324"),
            ("soil_thermal_cap", "1230000.0"),
            ("soil_thermal_cond", "0.22"),
            ("stdev_topog_index", "1.0"),
            ("water_table", "1.0"),
            ("wetness_under_soil", "1.0"),
        ]:
            self.add_setting(config, ["namelist:ideal_surface", name], value)
        """Add the albedo weights used with site forcing, the minimum LAI of
        snow-covered canopies and the soil thermal conductivity model, with
        the values previously used"""
        self.add_setting(
            config, ["namelist:jules_radiation", "wght_alb"], "0.0,0.5,0.0,0.5"
        )
        self.add_setting(
            config,
            ["namelist:jules_snow", "lai_alb_lim_sn"],
            "1.0,1.0,0.1,0.1,0.1",
        )
        self.add_setting(
            config, ["namelist:jules_soil", "soilhc_method"], "'johansen'"
        )
        """Move cansnowpft from jules_snow to cansnowpft_io in each instance
        of jules_pftparm, using the index of the instance's plant type"""
        cansnowpft = self.get_setting_value(
            config, ["namelist:jules_snow", "cansnowpft"]
        )
        if cansnowpft is not None:
            values = []
            for value in cansnowpft.replace("\n", "").split(","):
                value = value.strip().lstrip("=")
                if "*" in value:
                    count, value = value.split("*")
                    values += [value] * int(count)
                elif value:
                    values.append(value)
            for section in list(config.value.keys()):
                if not section.startswith("namelist:jules_pftparm("):
                    continue
                pft_name = self.get_setting_value(
                    config, [section, "pft_name_io"]
                )
                index = self.get_setting_value(
                    config,
                    ["namelist:jules_surface_types", pft_name.strip("'\"")],
                )
                if index is not None and 1 <= int(index) <= len(values):
                    value = values[int(index) - 1]
                else:
                    value = ".false."
                self.add_setting(config, [section, "cansnowpft_io"], value)
            self.remove_setting(config, ["namelist:jules_snow", "cansnowpft"])
        return config, self.reports
