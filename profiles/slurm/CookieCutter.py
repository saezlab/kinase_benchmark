#
# Based on lsf CookieCutter.py
#
import os
import json

d = os.path.dirname(__file__)
with open(os.path.join(d, "settings.json")) as fh:
    settings = json.load(fh)


class CookieCutter:

    SBATCH_DEFAULTS = settings['SBATCH_DEFAULTS']
    CLUSTER_NAME = settings['CLUSTER_NAME']
    CLUSTER_CONFIG = settings['CLUSTER_CONFIG']

    @staticmethod
    def get_cluster_option() -> str:
        cluster = CookieCutter.CLUSTER_NAME
        if cluster != "":
            return f"--cluster={cluster}"
        return ""