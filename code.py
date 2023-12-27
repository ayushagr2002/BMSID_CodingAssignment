import numpy as np
from icecream import ic

def rotate_point(point, angle, line_direction_cosines):
    # convert angle in degree to radians
    angle = np.deg2rad(angle)
    rotation_matrix = [
        [
            np.cos(angle) + line_direction_cosines[0] ** 2 * (1 - np.cos(angle)),
            line_direction_cosines[0] * line_direction_cosines[1] * (1 - np.cos(angle)) - line_direction_cosines[2] * np.sin(angle),
            line_direction_cosines[0] * line_direction_cosines[2] * (1 - np.cos(angle)) + line_direction_cosines[1] * np.sin(angle)
        ],
        [
            line_direction_cosines[0] * line_direction_cosines[1] * (1 - np.cos(angle)) + line_direction_cosines[2] * np.sin(angle),
            np.cos(angle) + line_direction_cosines[1] ** 2 * (1 - np.cos(angle)),
            line_direction_cosines[1] * line_direction_cosines[2] * (1 - np.cos(angle)) - line_direction_cosines[0] * np.sin(angle)
        ],
        [
            line_direction_cosines[0] * line_direction_cosines[2] * (1 - np.cos(angle)) - line_direction_cosines[1] * np.sin(angle),
            line_direction_cosines[1] * line_direction_cosines[2] * (1 - np.cos(angle)) + line_direction_cosines[0] * np.sin(angle),
            np.cos(angle) + line_direction_cosines[2] ** 2 * (1 - np.cos(angle))
        ]
    ]

    return np.matmul(np.array(rotation_matrix), np.array(point))

def get_fourth_atom_coordinates(first_atom_coordinates, second_atom_coordinates, third_atom_coordinates, bond_length, bond_angle, dihedral_angle):
    first_atom_coordinates = np.array(first_atom_coordinates)
    second_atom_coordinates = np.array(second_atom_coordinates)
    third_atom_coordinates = np.array(third_atom_coordinates)

    p = (second_atom_coordinates - first_atom_coordinates) 
    q = (third_atom_coordinates - second_atom_coordinates) 

    ABC_plane_normal = np.cross(p, q)
    ABC_plane_normal = ABC_plane_normal / np.linalg.norm(ABC_plane_normal)
    ic(ABC_plane_normal)

    unit_vector_along_q = q / np.linalg.norm(q)

    ic(unit_vector_along_q)

    initial_pos = unit_vector_along_q * bond_length

    rotated_pos = rotate_point(initial_pos, 180-bond_angle, ABC_plane_normal)

    final_pos = rotate_point(rotated_pos, dihedral_angle, unit_vector_along_q)

    # return final_pos
    return final_pos + third_atom_coordinates


if __name__ == '__main__':
    atom_coordinates = [[], [], []]
    atom_coordinates[0] = [0, 0, 0]
    atom_coordinates[1] = [1.453, 0, 0]
    atom_coordinates[2] = [1.959, -0.628, -1.30]
    bond_length = 1.3
    bond_angle = 120
    dihedral_angle = 30
    # for i in range(3):
    #     print("Enter Coordinates of atom {}:".format(i+1))
    #     x, y, z = map(float, input().split())
    #     atom_coordinates.append([x, y, z])

    # bond_length = float(input("Enter bond length: "))
    # bond_angle = float(input("Enter bond angle: "))
    # dihedral_angle = float(input("Enter dihedral angle: "))

    fourth_atom_coordinates = get_fourth_atom_coordinates(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2], bond_length, bond_angle, dihedral_angle)

    print("Coordinates of fourth atom are: {}".format(fourth_atom_coordinates))
    
