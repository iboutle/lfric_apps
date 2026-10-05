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

class vn32_t258(MacroUpgrade):
    # Upgrade macro for Issue #258 by Mike Whitall

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t258"

    def upgrade(self, config, meta_config=None):
        # Add settings

        # Blank macro needed due to meta-data changes.

        return config, self.reports


class vn32_t241(MacroUpgrade):
    # Upgrade macro for Issue#241 by Mike Whitall

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t241"

    def upgrade(self, config, meta_config=None):
        # Add settings

        # Microphysics settings
        nml = "namelist:microphysics"
        self.add_setting(config, [nml,"fix_casim_tidy"], ".false.")

        return config, self.reports


class vn32_t014(MacroUpgrade):
    # Upgrade macro for casim Issue#14 by Mike Whitall

    BEFORE_TAG = "vn3.2_t241"
    AFTER_TAG = "vn3.2_t014"

    def upgrade(self, config, meta_config=None):
        # Add settings

        nml = "namelist:microphysics"
        self.add_setting(config, [nml, "casim_inhom_rain"], ".false.")

        return config, self.reports
