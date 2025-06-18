from cmlibs.zinc.context import Context
from cmlibs.zinc.result import RESULT_OK
from ssvtools.query_structure import get_vagus_structure_maps, evaluate_branch_start_coordinates
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
        fieldmodule = region.getFieldmodule()
        structure_map, common_branch_map = get_vagus_structure_maps(fieldmodule)
        expected_structure_map = {
            'left vagus nerve': [None, ['left A thoracic cardiopulmonary branch of vagus nerve',
                                        'left B thoracic cardiopulmonary branch of vagus nerve',
                                        'left superior laryngeal nerve']],
            'left A branch of superior laryngeal nerve': ['left superior laryngeal nerve', []],
            'left A thoracic cardiopulmonary branch of vagus nerve': ['left vagus nerve', []],
            'left B thoracic cardiopulmonary branch of vagus nerve': ['left vagus nerve', []],
            'left superior laryngeal nerve': ['left vagus nerve', ['left A branch of superior laryngeal nerve']]}
        self.assertEqual(expected_structure_map, structure_map)
        expected_common_branch_map = {
            'left branch of superior laryngeal nerve': [
                'left A branch of superior laryngeal nerve'],
            'left thoracic cardiopulmonary branch of vagus nerve': [
                'left A thoracic cardiopulmonary branch of vagus nerve',
                'left B thoracic cardiopulmonary branch of vagus nerve']}
        self.assertEqual(common_branch_map, expected_common_branch_map)
        vagus_coordinates = fieldmodule.findFieldByName('vagus coordinates')
        branch_start_coordinates = evaluate_branch_start_coordinates(
            vagus_coordinates, structure_map['left vagus nerve'][1])
        expected_branch_start_coordinates = [
            ('left A thoracic cardiopulmonary branch of vagus nerve',
                [-0.0003328408968624165, 1.7493737182796838e-05, 0.3813658363063589],
                [-0.9106917613635638, 0.38295521752225525, 0.15487355215469165]),
            ('left B thoracic cardiopulmonary branch of vagus nerve',
                [0.0010020650864467912, -0.001554513928071552, 0.40674113258208205],
                [0.5645932995214857, -0.630854592172089, 0.5322188362608271]),
            ('left superior laryngeal nerve',
                [0.0008767095326078595, 0.00045882124897275354, 0.1315331326187551],
                [0.6662449900204912, 0.64287539363479, -0.3779270320200837])]
        self.assertEqual(len(branch_start_coordinates), len(expected_branch_start_coordinates))
        TOL = 1.0E-8
        for i in range(len(branch_start_coordinates)):
            expected_branch_name, expected_start_location, expected_direction = expected_branch_start_coordinates[i]
            branch_name, start_location, direction = branch_start_coordinates[i]
            self.assertEqual(branch_name, expected_branch_name)
            assertAlmostEqualList(self, start_location, expected_start_location, TOL)
            assertAlmostEqualList(self, direction, expected_direction, TOL)


if __name__ == "__main__":
    unittest.main()
