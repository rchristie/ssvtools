from cmlibs.zinc.context import Context
from cmlibs.zinc.result import RESULT_OK
from ssvtools.query_structure import query_branches
from testutils import assertAlmostEqualList
import os
import unittest


here = os.path.abspath(os.path.dirname(__file__))


class SSVToolsTestCase(unittest.TestCase):

    def test_query_branches(self):
        """
        Get positions and orientations of branches from a test subject-specific-vagus SSV.
        Uses a simple test scaffold, but any vagus scaffold output by SPARC Mapping Tools
        Scaffold Creator, including from REVA data should work.
        """
        data_file = os.path.join(here, "resources", "vagus_test_scaffold1.exf")
        context = Context("Test")
        region = context.getDefaultRegion()
        self.assertEqual(RESULT_OK, region.readFile(data_file))
        branch_info, common_branch_map = query_branches(region)
        expected_branch_info = {
            'left A thoracic cardiopulmonary branch of vagus nerve': (
                [-0.0003328408968624165, 1.7493737182796838e-05, 0.3813658363063589],
                [-0.9097938227559508, 0.38257762513656673, 0.16096447067591738]),
            'left B thoracic cardiopulmonary branch of vagus nerve': (
                [0.0010020650864467912, -0.001554513928071552, 0.40674113258208205],
                [0.5530901132994802, -0.6180013793215057, 0.5587178373088555]),
            'left superior laryngeal nerve': (
                [0.0008767095326078595, 0.00045882124897275354, 0.1315331326187551],
                [0.6981983439472732, 0.6737079331530234, -0.24215014622330203])}
        self.assertEqual(len(branch_info), len(expected_branch_info))
        TOL = 1.0E-8
        for branch_name in branch_info.keys():
            start_location, start_direction = branch_info[branch_name]
            expected_start_location, expected_start_direction = expected_branch_info[branch_name]
            assertAlmostEqualList(self, start_location, expected_start_location, TOL)
            assertAlmostEqualList(self, start_direction, expected_start_direction, TOL)

        expected_common_branch_map = {
            'left branch of superior laryngeal nerve': [
                'left A branch of superior laryngeal nerve'],
            'left thoracic cardiopulmonary branch of vagus nerve': [
                'left A thoracic cardiopulmonary branch of vagus nerve',
                'left B thoracic cardiopulmonary branch of vagus nerve']}
        self.assertEqual(common_branch_map, expected_common_branch_map)


if __name__ == "__main__":
    unittest.main()
