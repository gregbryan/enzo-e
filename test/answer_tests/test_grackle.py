import os
import numpy as np
import yt

from answer_testing import \
    EnzoETest, \
    ytdataset_test, \
    assert_array_rel_equal, \
    uses_grackle

_base_file = os.path.basename(__file__)

# Set test tolerance based on compile precision
use_double = os.environ.get("USE_DOUBLE", "false").lower() == "true"
yt.mylog.info(f"{_base_file}: use_double = {use_double}")
if use_double:
    decimals = 12
else:
    decimals = 6


def _get_decimals(field_name):
    # returns infinity (require an exact match), when we are looking at the
    # location of a minimum/maximum.
    assert isinstance(field_name, tuple)
    if field_name[1].split(':')[-1] not in ["min", "max", "mean", "std_dev"]:
        return np.inf
    elif use_double:
        return 14
    else:
        return 6

@uses_grackle
class TestGrackleGeneral(EnzoETest):
    parameter_file = "Grackle/method_grackle_general.in"
    max_runtime = 60
    ncpus = 4

    @ytdataset_test(assert_array_rel_equal, decimals = _get_decimals)
    def test_grackle_general(self):
        """
        Compares a variety of summary statistics.

        This was adapted from an earlier test.
        """
        ds = yt.load("GeneralGrackle-500.00/GeneralGrackle-500.00.block_list")
        ad = ds.all_data()

        quan_entry_sets = [
            ('min_location', ('min', 'min_xloc', 'min_yloc', 'min_zloc'), {}),
            ('max_location', ('max', 'max_xloc', 'max_yloc', 'max_zloc'), {}),
            ('weighted_standard_deviation',
             ('std_dev', 'mean'), {'weight' : ("gas", "cell_volume")}),
        ]

        data = {}
        for field in ds.field_list:
            for derived_quantity, rslt_names, kwargs in quan_entry_sets:
                quan_func = getattr(ad.quantities,derived_quantity)
                rslts = quan_func(field, **kwargs)
                num_rslts = len(rslts)
                assert num_rslts == len(rslt_names)
                for i in range(num_rslts):
                    data[f'{field[1]}:{rslt_names[i]}'] = rslts[i]
        return data
