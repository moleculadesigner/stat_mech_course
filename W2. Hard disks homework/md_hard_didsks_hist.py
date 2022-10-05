from copy import deepcopy
import random
from sys import argv
import math
import matplotlib.pyplot as plt
import numpy as np


def vec_diff(vec_a, vec_b):
    return map(
        lambda p: p[1] - p[0],
        zip(vec_a, vec_b)
    )

def vec_sq_distance(vec_a, vec_b):
    return sum(
        map(
            lambda x: pow(x, 2),
            vec_diff(vec_a, vec_b)
        )
    )

def dot_product(vec_a, vec_b):
    return sum(
        map(
            lambda p: p[0]*p[1],
            zip(vec_a, vec_b)
        )
    )


def wall_time(pos_a, vel_a, sigma):
    if vel_a > 0.0:
        del_t = (1.0 - sigma - pos_a) / vel_a
    elif vel_a < 0.0:
        del_t = (pos_a - sigma) / abs(vel_a)
    else:
        del_t = float('inf')
    return del_t


def pair_time(pos_a, vel_a, pos_b, vel_b, sigma):
    del_x = vec_diff(pos_a, pos_b)
    del_x_sq = vec_sq_distance(pos_a, pos_b)
    
    del_v = vec_diff(vel_a, vel_b)
    del_v_sq = vec_sq_distance(vel_a, vel_b)
    
    scal = dot_product(del_v, del_x)
    
    Upsilon = pow(scal, 2) - del_v_sq*(del_x_sq - 4.0*pow(sigma, 2))
    if Upsilon > 0.0 and scal < 0.0:
        del_t = -(scal + math.sqrt(Upsilon))/del_v_sq
    else:
        del_t = float('inf')
    return del_t


def run_multi_events(
        n_events,
        positions,
        velocities,
        sigma=0.1197,
    ):

    t = .0
    samples = []
    singles = [
        (particle, axis)
        for particle in range(4)
        for axis in range(2)
    ]
    pairs = [
        (p1, p2)
        for p1 in range(4)
        for p2 in range(4)
        if p2 > p1
    ]
    
    for event in range(n_events):
        if event % 100_000 == 0:
            print(event)
            print(*map(lambda v: f"[{v[0]:+.2f} {v[1]:+.2f}]", positions))
            print(*map(lambda v: f"[{v[0]:+.2f} {v[1]:+.2f}]", velocities))
    
        wall_times = [
            wall_time(
                positions[particle][axis],
                velocities[particle][axis],
                sigma
            )
            for particle, axis  in singles
        ]
        pair_times = [
            pair_time(
                positions[p1],
                velocities[p1],
                positions[p2],
                velocities[p2],
                sigma
            )
            for p1, p2 in pairs
        ]
        next_event = min(wall_times + pair_times)
    
        t_previous = t
        for inter_times in range(int(t + 1), int(t + next_event + 1)):
            del_t = inter_times - t_previous
            for particle, axis in singles:
                positions[particle][axis] += velocities[particle][axis]*del_t
            t_previous = inter_times
            
            samples.append(positions[0][0])

        t += next_event
        del_t = t - t_previous
        # why?
        for particle, axis in singles:
            positions[particle][axis] += velocities[particle][axis]*del_t
        if min(wall_times) < min(pair_times):
            collision_disk, direction = singles[wall_times.index(next_event)]
            velocities[collision_disk][direction] *= -1.0
        else:
            a, b = pairs[pair_times.index(next_event)]
            del_x = [positions[b][0] - positions[a][0], positions[b][1] - positions[a][1]]
            abs_x = math.sqrt(del_x[0] ** 2 + del_x[1] ** 2)
            e_perp = [c / abs_x for c in del_x]
            del_v = [velocities[b][0] - velocities[a][0], velocities[b][1] - velocities[a][1]]
            scal = dot_product(del_v, e_perp)
            for k in range(2):
                velocities[a][k] += e_perp[k] * scal
                velocities[b][k] -= e_perp[k] * scal

    return samples

pos = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]
vel = [[0.21, 0.12], [0.71, 0.18], [-0.23, -0.79], [0.78, 0.1177]]
n_events = int(argv[1])
samples = np.asarray(run_multi_events(n_events, pos, vel))

fig, ax = plt.subplots(figsize=(8,7))
ax.hist(
    samples,
    bins=100,
    color="k",
)
ax.set_xlabel('x')
ax.set_ylabel('frequency')
ax.set_title('MD sampling: x coordinate histogram (density eta=0.18)')
fig.savefig('md_disks_histo.png')
