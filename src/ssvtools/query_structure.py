"""
Utility functions for querying structure of SPARC subject-specific nerve scaffolds,
including trunk information, branch positions and orientations, and level marker locations.
"""
from cmlibs.maths.vectorops import normalize
from cmlibs.utils.zinc.field import get_group_list
from cmlibs.zinc.field import Field
from cmlibs.zinc.region import Region
import logging
import re


logger = logging.getLogger(__name__)


# following keyword is present in names of all annotations marking part of the trunk
trunk_keyword = 'vagus nerve'
# exception
non_trunk_keyword = 'of vagus nerve'
# following keywords are part of names of all annotations marking branches (without trunk_keywords)
branch_keywords = ['branch', 'nerve']


def get_trunk_group(fieldmodule):
    """
    Get the Zinc Group containing the vagus trunk part of the model.
    :param fieldmodule: Fieldmodule for a Zinc Region
    :return: Zinc Group for left/right vagus nerve or None if none.
    """
    for side in ('left', 'right'):
        group_name = side + ' ' + trunk_keyword
        # a group contains a set of model objects (elements, nodes), but is also
        # a field returning True at any location inside those parts, False outside.
        group = fieldmodule.findFieldByName(group_name).castGroup()
        if group.isValid():
            return group
    return None


def query_branches(region: Region):
    """
    Get information about position and orientation of branches off the vagus trunk.
    :param region: Zinc Region containing SPARC subject-specific vagus scaffold (EX file read into region).
    :return: dict mapping from branch name to start (vagus coordinates, body direction left, anterior, down)
    excluding common branches, dict mapping common branch names to list of variant branch names within them.
    """
    fieldmodule = region.getFieldmodule()
    # coordinates = fieldmodule.findFieldByName('coordinates').castFiniteElement()
    straight_coordinates = fieldmodule.findFieldByName('straight coordinates').castFiniteElement()
    vagus_coordinates = fieldmodule.findFieldByName('vagus coordinates').castFiniteElement()
    mesh3d = fieldmodule.findMeshByDimension(3)
    derivative_xi1 = mesh3d.getChartDifferentialoperator(1, 1)
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    trunk_group = get_trunk_group(fieldmodule)
    trunk_nodeset_group = trunk_group.getNodesetGroup(nodes)
    if not trunk_group:
        logger.error('query_branches: Missing vagus trunk group')
        return {}
    # all branches and other annotations are groups within the region
    # determine what is a branch by name matching.
    branch_info = {}
    group_list = get_group_list(fieldmodule)
    # vagus scaffolds have common groups for all branches with variant A/B/C names
    # first determine these groups; a map from common group name to list of variant group names
    common_group_map = {}
    for group in group_list:
        group_name = group.getName()
        common_name = re.sub(r'\b[A-Z]\b\s?', '', group_name).strip()
        if common_name != group_name:
            variant_group_list = common_group_map.get(common_name)
            if not variant_group_list:
                variant_group_list = []
            variant_group_list.append(group_name)
            common_group_map[common_name] = variant_group_list
    fieldcache = fieldmodule.createFieldcache()
    for group in group_list:
        group_name = group.getName()
        if common_group_map.get(group_name):
            continue
        if (trunk_keyword in group_name) and not (non_trunk_keyword in group_name):
            continue
        # branches have 3-D elements in them
        mesh_group = group.getMeshGroup(mesh3d)
        if (not mesh_group.isValid()) or (mesh_group.getSize() == 0):
            continue
        if any(branch_keyword in group_name for branch_keyword in branch_keywords):
            element = mesh_group.createElementiterator().next()
            eft = element.getElementfieldtemplate(vagus_coordinates, -1)
            node = element.getNode(eft, 1)
            if not trunk_nodeset_group.containsNode(node):
                continue  # parent is not the trunk
            # this is the centre of the branch at its root
            fieldcache.setMeshLocation(element, [0.0, 0.5, 0.5])
            # Note: 3rd component of material coordinate is unit coordinate down vagus trunk
            _, branch_vagus_coordinate = vagus_coordinates.evaluateReal(fieldcache, 3)
            _, d1 = straight_coordinates.evaluateDerivative(derivative_xi1, fieldcache, 3)
            # branch direction axes are (left, anterior, down the body)
            branch_direction = normalize(d1)
            branch_info[group_name] = (branch_vagus_coordinate, branch_direction)

    return branch_info, common_group_map
